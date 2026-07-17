"""Euler-FD Gamma-parts on the GATE systems (H2O, HF) -> per-pair-side tau values at
the exact gate geometries. H2O's H-part mixes O-H and H..H pairs; tau_HH(2.861) comes
from the H2+ curve tail (~-0.04, small) and tau_HO solves from the rest. The O-part
gives tau_OH directly. Same for HF (tau_HF, tau_FH). Onsite subtracted by the decoded
closed form at the system's actual Mulliken shell charges."""
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
GAM = {1: 0.8141717488, 8: 0.1454, 9: 0.567856801}
U = {1: [1.0062], 8: [1.4936, 0.8473], 9: [1.5589, 0.7561]}
NSH = {1: 1, 8: 2, 9: 2}
SYM = {1: "H", 8: "O", 9: "F"}
DELTA = 0.05

ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
H2O_XYZ = [[0.0, 0.0, 0.0],
           [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
           [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
CASES = [("H2O", [8, 1, 1], H2O_XYZ, 0),
         ("HF", [1, 9], [[0, 0, 0], [0, 0, 1.733]], 0)]


def es23(atoms, charge):
    r = oracle.run(atoms, charge=charge)
    return float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)", r["raw"], re.M).group(1))


def shell_q(zs, xyz, charge):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge)
    P, S = rec["state"]["P"], rec["S"]
    m = np.diag((P / 2.0) @ S)
    meta = []
    for at, z in enumerate(zs):
        meta.append((at, 0))
        if NSH[z] > 1:
            meta += [(at, 1)] * 3
    q = {}
    for mm, (at, l) in zip(m, meta):
        q[(at, l)] = q.get((at, l), 0.0) + 2 * mm
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


def pair_coef(q, zs, Rab, z_sel, pair):
    """coefficient of tau for Gamma_(z_sel)-carrying offsite terms on the given
    unordered atom pair: (1/6) * sum q_lA q_lB qat_A Gamma_A over ordered (A,B) in the
    pair with zs[A] == z_sel."""
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
            c += (1 / 6) * qa * qb * qat[at] * GAM[zs[at]]
    return c


P0 = params.parse()
out = {}
try:
    for name, zs, xyz, chg in CASES:
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        q = shell_q(zs, xyz, chg)
        nat = len(zs)
        Rab = [[float(np.linalg.norm(np.array(xyz[a]) - np.array(xyz[b])))
                for b in range(nat)] for a in range(nat)]
        base = es23(atoms, chg)
        for z_sel in sorted(set(zs)):
            g0 = GAM[z_sel]
            params.write_perturbed({(z_sel, 1, 6): g0 + DELTA})
            ep = es23(atoms, chg)
            params.write_perturbed({(z_sel, 1, 6): g0 - DELTA})
            em = es23(atoms, chg)
            params.write_perturbed({})
            part = g0 * (ep - em) / (2 * DELTA)
            off = part - e3_onsite_part(q, zs, z_sel)
            out[(name, z_sel)] = (q, zs, Rab, off)
            print(f"  {name} G_{SYM[z_sel]}: E3-part {part:+.6f}  offsite {off:+.6f}",
                  flush=True)
finally:
    params.write_perturbed({})

# ---- solve taus at the gate geometries
print("\nper-pair-side tau at gate geometries:")
# HF: single pair (0,1): tau_H (Gamma_H side) and tau_F
for z_sel, lab in ((1, "tau_H(HF,1.733)"), (9, "tau_F(HF,1.733)")):
    q, zs, Rab, off = out[("HF", z_sel)]
    c = pair_coef(q, zs, Rab, z_sel, (0, 1))
    print(f"  {lab} = {off / c:+.5f}   (coef {c:+.6f})")

# H2O O-part: O's partners are the two H's (equivalent O-H pairs)
q, zs, Rab, off = out[("H2O", 8)]
c = pair_coef(q, zs, Rab, 8, (0, 1)) + pair_coef(q, zs, Rab, 8, (0, 2))
tau_OH = off / c
print(f"  tau_O(H2O O-H,1.809) = {tau_OH:+.5f}   (coef {c:+.6f})")

# H2O H-part: two O-H pairs (unknown tau_HO) + one H..H pair (tau from H2+ tail)
q, zs, Rab, off = out[("H2O", 1)]
c_oh = pair_coef(q, zs, Rab, 1, (0, 1)) + pair_coef(q, zs, Rab, 1, (0, 2))
c_hh = pair_coef(q, zs, Rab, 1, (1, 2))
TAU_HH_2861 = -0.040   # H2+ curve tail interpolation (R = 2.861), small; labeled
tau_HO = (off - c_hh * TAU_HH_2861) / c_oh
print(f"  tau_H(H2O O-H,1.809) = {tau_HO:+.5f}   (coef {c_oh:+.6f}; H..H via H2+ tail "
      f"{TAU_HH_2861} x {c_hh:+.6f})")
print(f"\n  R(O-H) = {Rab[0][1]:.4f}  R(H..H) = {Rab[1][2]:.4f}")
