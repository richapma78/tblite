"""eeqbc.py -- OUR implementation of the EEQ(BC) charge model, gated against the v1 oracle.

The model (Froitzheim, Mueller, Hansen, Grimme, J. Chem. Phys. 2025, 162, 214109) as implemented
in multicharge's eeqbc2025 -- which the g-xTB v1 parameter file matches BYTE-FOR-BYTE across all
103 elements x 8 columns (extract_tables.py proves it on every run of the pipeline). Chain:

  CN_i     = sum_j 0.5*(1+erf(-kcn*(r-rc)/rc^norm_exp)),   rc = rcov_i + rcov_j   (erf counting)
  qloc_i   = sum_j count_ij * (EN_j - EN_i) + Q/nat        (EN-weighted CN + spread charge)
  c_ij     = sqrt(cap_i*cap_j) * 0.5*(1+erf(-kbc*(r-rvdw)/rvdw))          (bond capacitance)
  C        = Maxwell capacitance matrix (diag = sum of pair caps), bordered
  rad_i    = rad_i * (1 - kcnrad*CN_i/avg_cn_i^norm_exp)                  (CN-scaled width)
  A_ij     = erf(r/sqrt(radi^2+radj^2))/r * C_ij ;  A_ii = (eta_i + kqeta_i*qloc_i
             + sqrt(2/pi)/rad_i)*C_ii + 1 ;  bordered with the charge constraint
  x_i      = -chi_i + kcnchi_i*CN_i + kqchi_i*qloc_i ;  x_{n+1} = Q ;  RHS = C @ x
  solve A q = RHS  ->  partial charges q

Constants from multicharge's new_eeqbc2025_model: kcn=2.0, norm_exp=0.75, kbc=0.60, kcnrad=0.14,
cutoff=25 Bohr. EN = Pauling/3.98 with the actinide patches. rvdw = the pair table's ANGSTROM
literal used against r in BOHR -- multicharge converts to Bohr at declaration and back to
Angstrom at use (get_vdw_rad(...)*autoaa); we reproduce the code as it runs, not as one might
wish it were. avg_cn is not in the v1 file; the published array is used (same vintage).

GATE (declared): reproduce the oracle's printed EEQ block -- CN, q_CN (the local charge) and q --
to <= 5e-4 (its print precision) on water AND a diverse probe set. --test runs it.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as _K

KCN, NORM_EXP, KBC, KCNRAD, CUTOFF = (_K.EEQBC["kcn"], _K.EEQBC["norm_exp"],
                                      _K.EEQBC["kbc"], _K.EEQBC["kcnrad"],
                                      _K.EEQBC["cutoff_bohr"])
SQRT2PI = math.sqrt(2.0 / math.pi)

_T = json.load(open(os.path.join(HERE, "data", "mctc-tables.json")))
_P = _T["eeqbc2025"]
_VDW = _T["vdwrad_packed_lower_angstrom"]


def _en_table():
    en = list(_T["pauling_en"])
    for z in range(90, 104):
        en[z - 1] = 1.30
    en[87 - 1] = 0.80
    en[89 - 1] = 1.00
    for z in (90, 91, 92, 95):
        en[z - 1] = 1.10
    for z in (93, 94, 97, 103):
        en[z - 1] = 1.20
    return [x / 3.98 for x in en]


_EN = _en_table()


def _rvdw(zi, zj):
    lo, hi = min(zi, zj), max(zi, zj)
    return _VDW[lo - 1 + hi * (hi - 1) // 2]


def charges(zs, xyz_bohr, charge=0):
    """zs: atomic numbers; xyz_bohr: (n,3) Bohr. Returns dict(cn, qloc, q)."""
    zs = list(zs)
    xyz = np.asarray(xyz_bohr, float)
    n = len(zs)
    p = {k: np.array([_P[k][z - 1] for z in zs]) for k in
         ("chi", "eta", "rad", "kcnchi", "kqchi", "kqeta", "cap", "cov_radii", "avg_cn")}
    rcov = 0.5 * p["cov_radii"]
    en = np.array([_EN[z - 1] for z in zs])

    r = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    np.fill_diagonal(r, np.inf)
    inside = r <= CUTOFF

    rc = rcov[:, None] + rcov[None, :]
    count = 0.5 * (1.0 + np.vectorize(math.erf)(-KCN * (r - rc) / rc**NORM_EXP)) * inside
    cn = count.sum(axis=1)
    qloc = (count * (en[None, :] - en[:, None])).sum(axis=1) + charge / n

    rvdw = np.array([[_rvdw(zi, zj) for zj in zs] for zi in zs], float)
    cap = np.sqrt(p["cap"][:, None] * p["cap"][None, :])
    cpair = cap * 0.5 * (1.0 + np.vectorize(math.erf)(-KBC * (r - rvdw) / rvdw)) * inside

    C = np.zeros((n + 1, n + 1))
    C[:n, :n] = -cpair
    np.fill_diagonal(C[:n, :n], cpair.sum(axis=1))
    C[n, n] = 1.0

    rad = p["rad"] * (1.0 - KCNRAD * cn / p["avg_cn"]**NORM_EXP)
    gam = 1.0 / np.sqrt(rad[:, None]**2 + rad[None, :]**2)
    kernel = np.vectorize(math.erf)(r * gam) / r

    A = np.zeros((n + 1, n + 1))
    A[:n, :n] = kernel * C[:n, :n]
    np.fill_diagonal(A[:n, :n],
                     (p["eta"] + p["kqeta"] * qloc + SQRT2PI / rad) * np.diag(C[:n, :n]) + 1.0)
    A[n, :], A[:, n] = 1.0, 1.0
    A[n, n] = 0.0

    x = np.empty(n + 1)
    x[:n] = -p["chi"] + p["kcnchi"] * cn + p["kqchi"] * qloc
    x[n] = charge
    q = np.linalg.solve(A, C @ x)[:n]
    return {"cn": cn, "qloc": qloc, "q": q}


def _gate(atoms, charge=0, uhf=0, tol=5e-4, label=""):
    """Compare against the oracle's printed EEQ block for one structure."""
    import oracle
    r = oracle.run(atoms, charge=charge, uhf=uhf)
    if not r["eeq"]:
        return None, f"{label}: oracle printed no EEQ block"
    BOHR = 1.8897261254578281
    zs = [_SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
    xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * BOHR
    ours = charges(zs, xyz, charge=charge)
    worst = {"cn": 0.0, "qloc": 0.0, "q": 0.0}
    for i, row in enumerate(r["eeq"]):
        worst["cn"] = max(worst["cn"], abs(ours["cn"][i] - row["cn"]))
        worst["qloc"] = max(worst["qloc"], abs(ours["qloc"][i] - row["q_cn"]))
        worst["q"] = max(worst["q"], abs(ours["q"][i] - row["q"]))
    ok = all(v <= tol for v in worst.values())
    return ok, (f"{label:14s} CN {worst['cn']:.1e}  qloc {worst['qloc']:.1e}  "
                f"q {worst['q']:.1e}  {'MATCH' if ok else 'MISMATCH'}")


_SYM2Z = {"h": 1, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "si": 14, "p": 15, "s": 16,
          "cl": 17, "br": 35, "i": 53, "pd": 46}

PROBES = {
    "water": ([("O", 0.0, 0.0, 0.117), ("H", 0.0, 0.757, -0.471), ("H", 0.0, -0.757, -0.471)], 0, 0),
    "HCl": ([("Cl", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.2746)], 0, 0),
    "CO": ([("C", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.128)], 0, 0),
    "NH4+": ([("N", 0.0, 0.0, 0.0), ("H", 1.0134, 0.0, 0.0), ("H", -0.3378, 0.9553, 0.0),
              ("H", -0.3378, -0.4776, -0.8273), ("H", -0.3378, -0.4776, 0.8273)], 1, 0),
    "AcCl": ([("C", -1.428, 0.257, 0.0), ("C", 0.062, 0.400, 0.0), ("O", 0.616, 1.472, 0.0),
              ("Cl", 0.930, -1.147, 0.0), ("H", -1.877, 1.249, 0.0), ("H", -1.699, -0.304, 0.894),
              ("H", -1.699, -0.304, -0.894)], 0, 0),
    "PdCl2": ([("Pd", 0.0, 0.0, 0.0), ("Cl", 2.31, 0.0, 0.0), ("Cl", -2.31, 0.0, 0.0)], 0, 0),
}


if __name__ == "__main__":
    fails = 0
    for label, (atoms, chg, uhf) in PROBES.items():
        ok, msg = _gate(atoms, charge=chg, uhf=uhf, label=label)
        print("  " + msg)
        fails += 0 if ok else 1
    if fails:
        print(f"\nEEQBC GATE: {fails}/{len(PROBES)} MISMATCH -- the model chain is not yet the oracle's")
        sys.exit(1)
    print("\nEEQBC GATE PASSED: our charges ARE the oracle's, across the probe set")
