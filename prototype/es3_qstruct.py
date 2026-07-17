"""The charge-state scan: the q-structure of the offsite ES3, decoded by falsification.
Five states of the O-H pair at fixed geometries -- OH- / OH(radical) / H2O / H2O+ /
H2O2+ -- each Euler-FD'd per element (Gamma_X * dES23/dGamma_X = X-part, exact),
onsite subtracted by the decoded closed form. MODELS for the Gamma_X-carrying offsite
term (per ordered pair A-B, X = element of A):
  M1: (1/6) q_lA q_lB qat_A tau G_A      (Eq-129b as extracted)
  M2: (1/6) q_lA q_lB qat_B tau G_A      (partner-charge weighting)
  M3: (1/6) q_lA q_lB (qat_A+qat_B)/2 tau G_A
PRE-DECLARED: a model WINS if one tau per pair-side fits all states with spread <= 5%;
otherwise bank per-state taus and the structural failure. delta_FD: H 0.05, O 0.02
(Gamma_O = 0.145 is small; keep the relative step sane)."""
import math
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

BOHR = K.BOHR
K3GS, K3GP = -0.05708, -0.06679
GAM = {1: 0.8141717488, 8: 0.1454}
U = {1: [1.0062], 8: [1.4936, 0.8473]}
NSH = {1: 1, 8: 2}
SYM = {1: "H", 8: "O"}
DELTA = {1: 0.05, 8: 0.02}
TAU_HH_2861 = -0.040   # H2+ curve tail, H..H in the water frame

ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
W = [[0.0, 0.0, 0.0],
     [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
     [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
CASES = [("OH-", [8, 1], [[0, 0, 0], [0, 0, 1.83]], -1, 0),
         ("OH.", [8, 1], [[0, 0, 0], [0, 0, 1.83]], 0, 1),
         ("H2O", [8, 1, 1], W, 0, 0),
         ("H2O+", [8, 1, 1], W, 1, 1),
         ("H2O2+", [8, 1, 1], W, 2, 0)]


def es23(atoms, charge, uhf):
    r = oracle.run(atoms, charge=charge, uhf=uhf)
    return float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)", r["raw"], re.M).group(1))


def shell_q(zs, xyz, charge, uhf):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge, uhf=uhf)
    st = rec["state"]
    P = st["P"] if "P" in st else st["P_a"] + st["P_b"]
    S = rec["S"]
    net = np.diag(P @ S)
    meta = []
    for at, z in enumerate(zs):
        meta.append((at, 0))
        if NSH[z] > 1:
            meta += [(at, 1)] * 3
    q = {}
    for nn, (at, l) in zip(net, meta):
        q[(at, l)] = q.get((at, l), 0.0) + nn
    for (at, l) in list(q):
        q[(at, l)] = K.REFOCC[zs[at]].get(l, 0.0) - q[(at, l)]
    return q


def e3_onsite_part(q, zs, z_sel):
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    e3 = 0.0
    for (at, la), qa in q.items():
        if zs[at] != z_sel:
            continue
        z = zs[at]
        for lb in range(NSH[z]):
            qb = q[(at, lb)]
            ta = -1.0 / (2 * U[z][la] ** 2)
            tb = -1.0 / (2 * U[z][lb] ** 2)
            ga = (K3GS if la == 0 else K3GP) * GAM[z]
            gb = (K3GS if lb == 0 else K3GP) * GAM[z]
            e3 += (1 / 6) * qa * qb * qat[at] * (ta * ga + tb * gb)
    return e3


def coef(q, zs, z_sel, pair, model):
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    c = 0.0
    for (at, la), qa in q.items():
        if zs[at] != z_sel or at not in pair:
            continue
        for (bt, lb), qb in q.items():
            if bt == at or bt not in pair:
                continue
            w = qat[at] if model == 1 else (qat[bt] if model == 2
                                            else 0.5 * (qat[at] + qat[bt]))
            c += (1 / 6) * qa * qb * w * GAM[zs[at]]
    return c


P0 = params.parse()
data = []
try:
    for name, zs, xyz, chg, uhf in CASES:
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        try:
            q = shell_q(zs, xyz, chg, uhf)
        except Exception as e:
            print(f"  {name}: restart FAILED ({type(e).__name__}: {e}) -- skipped",
                  flush=True)
            continue
        qat = {}
        for (at, l), v in q.items():
            qat[at] = qat.get(at, 0.0) + v
        parts = {}
        for z_sel in sorted(set(zs)):
            g0 = GAM[z_sel]
            dl = DELTA[z_sel]
            params.write_perturbed({(z_sel, 1, 6): g0 + dl})
            ep = es23(atoms, chg, uhf)
            params.write_perturbed({(z_sel, 1, 6): g0 - dl})
            em = es23(atoms, chg, uhf)
            params.write_perturbed({})
            part = g0 * (ep - em) / (2 * dl)
            parts[z_sel] = part - e3_onsite_part(q, zs, z_sel)
        data.append((name, zs, q, parts))
        print(f"  {name}: qat " + " ".join(f"{SYM[zs[a]]}{a}:{qat[a]:+.3f}"
                                           for a in sorted(qat))
              + "   offsite G_O-part {:+.6f}  G_H-part {:+.6f}".format(
                  parts[8], parts[1]), flush=True)
finally:
    params.write_perturbed({})

print("\nimplied tau per model (O-side and H-side of the O-H pair):")
for model in (1, 2, 3):
    print(f"  MODEL M{model}:")
    for side, z_sel in (("tau_O", 8), ("tau_H", 1)):
        vals = []
        for name, zs, q, parts in data:
            pairs_oh = [(0, a) for a in range(1, len(zs))]
            c = sum(coef(q, zs, z_sel, pr, model) for pr in pairs_oh)
            off = parts[z_sel]
            if z_sel == 1 and len(zs) == 3:      # subtract the H..H piece
                c_hh = coef(q, zs, 1, (1, 2), model)
                off = off - c_hh * TAU_HH_2861
            if abs(c) < 1e-7:
                vals.append((name, None))
                continue
            vals.append((name, off / c))
        got = [v for _n, v in vals if v is not None]
        spread = (max(got) - min(got)) / max(abs(np.mean(got)), 1e-9) if got else 0
        line = "  ".join(f"{n}:{v:+.3f}" if v is not None else f"{n}:--"
                         for n, v in vals)
        print(f"    {side}: {line}   spread {100 * spread:.0f}%")
