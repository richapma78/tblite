"""The fit-free ES3 onsite validation ladder.
Measured once on He+: k3Gs*Gamma_He = 6*E3(He+) = 0.6676. Everything after is
PREDICTION, no fitting:
  H-  (s-only anion, q=-1): E3 = -(k3Gs*Gamma_H)/6, Gamma_H = L1[6](H) from the file,
      k3Gs = 0.6676/Gamma_He(file). ES2 = (1/2)*s_H*U_H. Printed ES2+3 is predicted.
  Li- (s-anion): same s-channel prediction with Li's file constants.
  F-  (p-anion): tests the p-channel: E3 ~ -(1/6)*[k3Gp*Gamma_F]*qp^3 + s-p cross with
      the harmonic-mean tau; k3Gp candidate = G2[5] = 0.4667. Mulliken shell charges
      from the restart (not assumed integer).
Gates (pre-declared): H-/Li- printed ES2+3 within 0.01 of prediction = s-channel
CONFIRMED; F- within 0.02 decides whether G2[5] is k3Gamma_p."""
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

P0 = params.parse()


def L1(z, col):
    # L1 row = first row after the element header; parse() stores rows as "shells"
    # offset-shifted -- read the pristine directly to be layout-proof
    import os
    lines = open(os.path.join(os.path.dirname(params.__file__), "data",
                              "gxtb_parameters.pristine")).read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == [str(z)])
    return float(lines[idx[idx.index(hdr) + 1]].split()[col])


def L6(z, col=0):
    import os
    lines = open(os.path.join(os.path.dirname(params.__file__), "data",
                              "gxtb_parameters.pristine")).read().splitlines()
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == [str(z)])
    return float(lines[idx[idx.index(hdr) + 6]].split()[col])


G_H, G_HE, G_LI, G_F = L1(1, 6), L1(2, 6), L1(3, 6), L1(9, 6)
print(f"Gamma column L1[6]: H {G_H}  He {G_HE}  Li {G_LI}  F {G_F}")
K3GS = 0.6676 / G_HE
print(f"k3Gs = 6*E3(He+)/Gamma_He = {K3GS:.5f}")

SR = {1: 0.4726, 3: 0.08247 * 1 + 0.0920, 9: 0.08247 * 7 + 0.0920}


def es23(atoms, charge=0, uhf=0):
    r = oracle.run(atoms, charge=charge, uhf=uhf)
    m = re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)", r["raw"], re.M)
    return float(m.group(1))


# ---- H- and Li- (s-channel predictions)
for name, z, sym in (("H-", 1, "H"), ("Li-", 3, "Li")):
    U = L6(z)
    G = L1(z, 6)
    es2 = 0.5 * SR[z] * U
    e3 = -(K3GS * G) / 6.0
    pred = es2 + e3
    try:
        got = es23([(sym, 0, 0, 0)], charge=-1, uhf=0)
        print(f"{name:3s}: predicted ES2+3 {pred:+.5f} (ES2 {es2:+.5f} + E3 {e3:+.5f})  "
              f"printed {got:+.5f}   d {got - pred:+.5f}"
              f"   {'CONFIRMED' if abs(got - pred) <= 0.01 else 'MISS'}", flush=True)
    except Exception as e:
        print(f"{name}: run failed ({e})")

# ---- F- (p-channel test, real Mulliken shell charges)
rec = fock_recon.fock_ao([("F", 0, 0, 0)], charge=-1)
P = rec["state"]["P"]
S = rec["S"]
m = np.diag((P / 2.0) @ S)
qs = K.REFOCC[9][0] - 2 * m[0]
qp = K.REFOCC[9][1] - 2 * (m[1] + m[2] + m[3])
qA = qs + qp
Us, Up = L6(9, 0), L6(9, 1)
print(f"\nF-: q_s {qs:+.4f}  q_p {qp:+.4f}")
g_ss, g_pp = Us, Up
g_sp = 2 * Us * Up / (Us + Up)
es2 = 0.5 * SR[9] * (qs * qs * g_ss + qp * qp * g_pp + 2 * qs * qp * g_sp)


def tau_on(Ua, Ub):
    return 2 * Ub * Ub / (Ua + Ub) ** 2


K3GP = 0.4667  # G2[5] candidate
GS_F = K3GS * G_F
GP_F = K3GP * G_F
e3 = 0.0
for (qa, Ua, Ga) in ((qs, Us, GS_F), (qp, Up, GP_F)):
    for (qb, Ub, Gb) in ((qs, Us, GS_F), (qp, Up, GP_F)):
        e3 += (1.0 / 6.0) * qa * qb * (qA * tau_on(Ua, Ub) * Ga
                                       + qA * tau_on(Ub, Ua) * Gb)
pred = es2 + e3
got = float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)",
                      rec["state"]["raw"], re.M).group(1))
print(f"F-: predicted ES2+3 {pred:+.5f} (ES2 {es2:+.5f} + E3 {e3:+.5f})  "
      f"printed {got:+.5f}   d {got - pred:+.5f}"
      f"   {'p-channel CONFIRMED (G2[5]=k3Gp)' if abs(got - pred) <= 0.02 else 'MISS'}")
