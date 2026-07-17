"""The final offsite decomposition: residue = measured_offsite_part - fold(c=1), then
tau_Q = residue / Q-weighted coefficient. PRE-DECLARED: the offsite ES3 is DECODED if
tau_Q is state-invariant per pair-side (<= 10% where signals are above floor) and the
H2+ -> H3+ transfer closes (the pass-76 20% miss should dissolve if it was fold
contamination). All offsite parts below are MEASURED (passes 76-78); only charges are
recomputed here."""
import math
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402

BOHR = K.BOHR
GAM = {1: 0.8141717488, 2: 5.0001493881, 8: 0.1454}
U = {1: [1.0062], 2: [0.6537], 8: [1.4936, 0.8473]}
NSH = {1: 1, 2: 1, 8: 2}
SYM = {1: "H", 2: "He", 8: "O"}

OFF = {  # measured offsite Gamma-parts (passes 76-78)
    ("H2+@1.2", 1): +0.074258, ("H2+@1.5", 1): +0.044828, ("H2+@1.8", 1): +0.021621,
    ("H2+@2.1", 1): +0.007286, ("H2+@2.5", 1): -0.000166, ("H2+@3.0", 1): -0.001797,
    ("H3+@1.5", 1): +0.04864, ("H3+@1.65", 1): +0.03547, ("H3+@1.9", 1): +0.01820,
    ("HeH+@1.46", 1): +0.02900, ("HeH+@1.46", 2): +0.01274,
    ("HeH+@2.0", 1): -0.00070, ("HeH+@2.0", 2): +0.00232,
    ("HeH+@2.5", 1): -0.00531, ("HeH+@2.5", 2): -0.00005,
    ("OH-@1.83", 1): -0.00527, ("OH-@1.83", 8): +0.01258,
    ("H2O+", 1): +0.022837, ("H2O+", 8): +0.001400,
}


def tri(R):
    return [[0, 0, 0], [R, 0, 0], [R / 2, R * math.sqrt(3) / 2, 0]]


ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
W = [[0.0, 0.0, 0.0],
     [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
     [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
CASES = {}
for R in (1.2, 1.5, 1.8, 2.1, 2.5, 3.0):
    CASES[f"H2+@{R}"] = ([1, 1], [[0, 0, 0], [0, 0, R]], 1, 1)
for R in (1.5, 1.65, 1.9):
    CASES[f"H3+@{R}"] = ([1, 1, 1], tri(R), 1, 0)
for R in (1.46, 2.0, 2.5):
    CASES[f"HeH+@{R}"] = ([2, 1], [[0, 0, 0], [0, 0, R]], 1, 0)
CASES["OH-@1.83"] = ([8, 1], [[0, 0, 0], [0, 0, 1.83]], -1, 0)
CASES["H2O+"] = ([8, 1, 1], W, 1, 1)


def shell_q(zs, xyz, charge, uhf):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge, uhf=uhf)
    st = rec["state"]
    P = st["P"] if "P" in st else st["P_a"] + st["P_b"]
    net = np.diag(P @ rec["S"])
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


def fold_part(q, zs, Rab, z_sel):
    """c=1 offsite-KO fold: Gamma_zsel * sum_{at of z_sel} qat * dES2_off/dU_l."""
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    s = 0.0
    for (at, la), qa in q.items():
        if zs[at] != z_sel:
            continue
        Us = U[z_sel][la]
        for (bt, lb), qb in q.items():
            if bt == at:
                continue
            g = 1.0 / (Rab[at][bt] + 0.5 * (1.0 / Us + 1.0 / U[zs[bt]][lb]))
            s += qat[at] * qa * qb * g * g / (2 * Us * Us)
    return GAM[z_sel] * s


def qcoef(q, zs, Q, z_sel):
    """(1/6) sum_ordered-cross-pairs q_lA q_lB * Q * Gamma_A over A of z_sel."""
    c = 0.0
    for (at, la), qa in q.items():
        if zs[at] != z_sel:
            continue
        for (bt, lb), qb in q.items():
            if bt == at:
                continue
            c += (1 / 6) * qa * qb * Q * GAM[z_sel]
    return c


print(f"{'state':12s} {'side':4s} {'offsite':>9s} {'fold(c=1)':>9s} {'residue':>9s} "
      f"{'tau_Q':>8s}")
taus = {}
for (name, z_sel), off in sorted(OFF.items()):
    zs, xyz, chg, uhf = CASES[name]
    q = shell_q(zs, xyz, chg, uhf)
    nat = len(zs)
    Rab = [[float(np.linalg.norm(np.array(xyz[a]) - np.array(xyz[b])))
            for b in range(nat)] for a in range(nat)]
    f = fold_part(q, zs, Rab, z_sel)
    res = off - f
    c = qcoef(q, zs, chg, z_sel)
    tq = res / c if abs(c) > 1e-7 else float("nan")
    fam = name.split("@")[0]
    taus.setdefault((fam, z_sel), []).append((name, tq))
    print(f"{name:12s} {SYM[z_sel]:4s} {off:+9.5f} {f:+9.5f} {res:+9.5f} {tq:+8.3f}")

print("\nH2+ vs H3+ transfer at matched R (tau_Q, the pass-76 20% question):")
h2 = dict((n.split("@")[1], t) for n, t in taus[("H2+", 1)])
h3 = dict((n.split("@")[1], t) for n, t in taus[("H3+", 1)])
for R in ("1.5",):
    print(f"  R={R}: H2+ tau_Q {h2[R]:+.3f}  H3+ tau_Q {h3[R]:+.3f}  "
          f"ratio {h3[R] / h2[R]:.3f}")
print("  (interpolate H2+ 1.65/1.9 from the ladder for the other two H3+ points)")
