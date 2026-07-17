"""ES3 completion sweep. Onsite form CONFIRMED (tau = -1/(2U^2), k3Gs = -0.05708;
H- fit-free match 4%). This script:
  1. k3Gp: solved from F- (one unknown, one equation), then the value is reported.
  2. Offsite (k3, k3x): E3_off targets = printed - oursES2 - E3_onsite on the H3+ ladder
     (1.5/1.65/1.9) and HeH+ ladder (1.46/2.0/2.5); solve the 2 unknowns by grid+refine
     over the closed form tau_off = d/dU[k3*Ubar^2*R*exp(-k3x*Ubar^2*R)]
     = k3*Ubar*R*e^(-k3x*Ubar^2*R)*(1 - k3x*Ubar^2*R); compare with candidate globals
     G2[0] = -1.8218, G2[1] = +0.3300, G2[5] = +0.4667.
All charges are real Mulliken shell charges from the restart."""
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402

BOHR = K.BOHR
K3GS = -0.05708
SR = {1: 0.4726, 2: 0.9218, 8: 0.08247 * 6 + 0.0920, 9: 0.08247 * 7 + 0.0920}
_LINES = open(os.path.join(os.path.dirname(params.__file__), "data",
                           "gxtb_parameters.pristine")).read().splitlines()
_IDX = [i for i, l in enumerate(_LINES) if l.strip()]


def row(z, off):
    hdr = next(i for i in _IDX if _LINES[i].split() == [str(z)])
    return [float(x) for x in _LINES[_IDX[_IDX.index(hdr) + off]].split()]


GAM = {z: row(z, 1)[6] for z in (1, 2, 8, 9)}
U = {z: row(z, 6) for z in (1, 2, 8, 9)}
NSH = {1: 1, 2: 1, 8: 2, 9: 2}
print(f"Gamma: { {z: round(v, 4) for z, v in GAM.items()} }")
print(f"U_s/U_p: { {z: [round(x, 4) for x in U[z][:NSH[z]]] for z in U} }")


def shell_q(rec, zs):
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
        q[(at, l)] = K.REFOCC[zs[at]].get(l, 2.0 if l == 0 else 0.0) - q[(at, l)]
    return q


def es2_ours(q, zs, Rab):
    es2 = 0.0
    for (at, l), qa in q.items():
        z = zs[at]
        for l2 in range(NSH[z]):
            if U[z][l] == 0.0 or U[z][l2] == 0.0:
                continue
            g2 = SR[z] * 2 * U[z][l] * U[z][l2] / (U[z][l] + U[z][l2])
            es2 += 0.5 * qa * q[(at, l2)] * g2
    for (at, l), qa in q.items():
        for (bt, l2), qb in q.items():
            if bt <= at or U[zs[at]][l] == 0.0 or U[zs[bt]][l2] == 0.0:
                continue
            gko = 1.0 / (Rab[at][bt] + 0.5 * (1.0 / U[zs[at]][l] + 1.0 / U[zs[bt]][l2]))
            es2 += qa * qb * gko
    return es2


def gam_l(z, l, k3gp):
    return (K3GS if l == 0 else k3gp) * GAM[z]


def e3_onsite(q, zs, k3gp):
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    e3 = 0.0
    for (at, la), qa in q.items():
        z = zs[at]
        for lb in range(NSH[z]):
            if U[z][la] == 0.0 or U[z][lb] == 0.0:
                continue
            qb = q[(at, lb)]
            ta = -1.0 / (2 * U[z][la] ** 2)
            tb = -1.0 / (2 * U[z][lb] ** 2)
            e3 += (1 / 6) * qa * qb * qat[at] * (ta * gam_l(z, la, k3gp)
                                                + tb * gam_l(z, lb, k3gp))
    return e3


def e3_offsite(q, zs, Rab, k3, k3x, k3gp, var=1):
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    e3 = 0.0
    for (at, la), qa in q.items():
        for (bt, lb), qb in q.items():
            if at == bt or U[zs[at]][la] == 0.0 or U[zs[bt]][lb] == 0.0:
                continue
            R = Rab[at][bt]
            ub = 0.5 * (U[zs[at]][la] + U[zs[bt]][lb])
            if var in (1, 2):        # bare Gamma_A offsite
                G = GAM[zs[at]]
            else:                    # onsite-style k3G_l * Gamma_A
                G = gam_l(zs[at], la, k3gp)
            if var in (1, 3):        # exponent with ub^2 (as SI-extracted)
                w = ub * ub
            else:                    # exponent with ub (alternate reading)
                w = ub
            ex = math.exp(-k3x * w * R)
            tau = k3 * ub * R * ex * (1 - k3x * w * R)
            e3 += (1 / 6) * qa * qb * (qat[at] * tau * G)
            # the qB*tau*Gamma_lB mirror comes from the (bt,at) loop iteration
    return e3


def printed(rec):
    return float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)",
                           rec["state"]["raw"], re.M).group(1))


# ---- 1. F- : solve k3Gp
rec = fock_recon.fock_ao([("F", 0, 0, 0)], charge=-1)
qF = shell_q(rec, [9])
tgt = printed(rec) - es2_ours(qF, [9], [[0.0]])
# e3_onsite is linear in k3gp: solve
a0 = e3_onsite(qF, [9], 0.0)
a1 = e3_onsite(qF, [9], 1.0) - a0
k3gp = (tgt - a0) / a1
print(f"\nF-: q_s {qF[(0, 0)]:+.4f} q_p {qF[(0, 1)]:+.4f}  target E3 {tgt:+.5f}  "
      f"=> k3Gp = {k3gp:+.5f}   (G2[5] = +0.4667 for comparison)")

# ---- 2. offsite ladders
SYM = {1: "H", 2: "He"}


def tri(Rr):
    return [[0, 0, 0], [Rr, 0, 0], [Rr / 2, Rr * math.sqrt(3) / 2, 0]]


targets = []
for name, zs, xyz, chg in [
        ("H3+@1.5", [1, 1, 1], tri(1.5), 1),
        ("H3+@1.65", [1, 1, 1], tri(1.65), 1),
        ("H3+@1.9", [1, 1, 1], tri(1.9), 1),
        ("HeH+@1.46", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 1),
        ("HeH+@2.0", [2, 1], [[0, 0, 0], [0, 0, 2.0]], 1),
        ("HeH+@2.5", [2, 1], [[0, 0, 0], [0, 0, 2.5]], 1)]:
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=chg)
    q = shell_q(rec, zs)
    nat = len(zs)
    Rab = [[float(np.linalg.norm(np.array(xyz[a]) - np.array(xyz[b])))
            for b in range(nat)] for a in range(nat)]
    off = printed(rec) - es2_ours(q, zs, Rab) - e3_onsite(q, zs, k3gp)
    targets.append((name, q, zs, Rab, off))
    print(f"  {name}: E3_offsite target {off:+.5f}")

VN = {1: "bare-G, exp(ub^2 R)", 2: "bare-G, exp(ub R)",
      3: "k3Gl*G, exp(ub^2 R)", 4: "k3Gl*G, exp(ub R)"}
overall = None
for var in (1, 2, 3, 4):
    best = None
    for k3 in np.arange(-8.0, 8.01, 0.2):
        if abs(k3) < 1e-9:
            continue
        for k3x in np.arange(0.02, 3.01, 0.04):
            r2 = sum((e3_offsite(q, zs, Rab, k3, k3x, k3gp, var) - off) ** 2
                     for _, q, zs, Rab, off in targets)
            if best is None or r2 < best[0]:
                best = (r2, k3, k3x)
    r2, k3, k3x = best
    for _ in range(4):
        for k3t in np.arange(k3 - 0.2, k3 + 0.2, 0.02):
            for k3xt in np.arange(max(0.01, k3x - 0.06), k3x + 0.06, 0.005):
                r2t = sum((e3_offsite(q, zs, Rab, k3t, k3xt, k3gp, var) - off) ** 2
                          for _, q, zs, Rab, off in targets)
                if r2t < r2:
                    r2, k3, k3x = r2t, k3t, k3xt
    rms = math.sqrt(r2 / len(targets))
    print()
    print(f"variant {var} ({VN[var]}): k3 = {k3:+.4f}  k3x = {k3x:.4f}  rms {rms:.5f}")
    for name, q, zs, Rab, off in targets:
        pred = e3_offsite(q, zs, Rab, k3, k3x, k3gp, var)
        print(f"    {name}: target {off:+.5f}  pred {pred:+.5f}  d {pred - off:+.5f}")
    if overall is None or rms < overall[0]:
        overall = (rms, var, k3, k3x)
print()
print(f"BEST: variant {overall[1]} k3 {overall[2]:+.4f} k3x {overall[3]:.4f} "
      f"rms {overall[0]:.5f}")
print("  (candidate globals: G2[0] -1.8218, G2[1] +0.3300, G2[5] +0.4667)")
