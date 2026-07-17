"""f2_scf.py -- the F2 SCF from decoded parts: the first p-shell self-consistent gate.

Assembly (every convention decided by measurement):
  H0    forward: file amplitudes x per-shell metric S~(k_s=0.730, k_p=1.205, kb=0.0632);
        diag = -L2_l; Pi ~ 1 (F's k_shp tiny)
  ES1   -1/2 S o (v+v), v_l = mu_l (corr -> 1 at the gate distance; measured curves say
        <= 0.5% there)
  ES2   -1/2 S o (v+v), v_l = sum_l' gamma2_onsite(l,l') q_l' with the DECODED kernel
        gamma2 = s(F)*harmonic(U_l, U_l') and q from the CURRENT density (self-consistent);
        offsite gamma2 omitted at the gate distance (q-weighted, < 0.005 -- bounded)
  X     the five-manifold skeleton per-AO/same-shell: v_i = 2*c_x*L5_l(i)*m_i, F = -1/2
        S o (v+v); plus the p-sigma object (banked closed form, R-lookup, labeled)
  ACP   analytic
GATE (pre-declared): all 8 eigenvalues within 0.02 Eh of the printed set at R = 3.8
(the unfitted p short pieces are small there); r_e reported as diagnostic (expect
0.03-0.05 misses from the short pieces).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapt  # noqa: E402
import basisq  # noqa: E402
import constants as K  # noqa: E402
import f2_stretch  # noqa: E402
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

P_ = params.parse()
eF = P_["element"][9]
kW_s, kW_p = P_["globals"]["g1"][0], P_["globals"]["g1"][1]
kdiat_sg, kdiat_pi = eF["l8"][0], eF["l8"][1]
L2 = [eF["shells"][0][0], eF["shells"][0][1]]
MU = [eF["shells"][5][0], eF["shells"][5][1]]
U = [eF["shells"][4][0], eF["shells"][4][1]]
L5 = [eF["shells"][3][0], eF["shells"][3][1]]
S_F = 0.08247 * 7 + 0.0920
CX = S_F / 9.59
KS, KP, KB = 0.730, 1.205, 0.0632
REF = K.REFOCC[9]
INCLUDE_SHORT = True   # f2_shortfit sets False to re-extract cleanly
LOFAO = [0, 1, 1, 1, 0, 1, 1, 1]                  # shell of each AO
ATOFAO = [0, 0, 0, 0, 1, 1, 1, 1]

Bq = basisq.parse()
_sh = [(l, np.array([p[0] for p in prims]), np.array([p[1] for p in prims]),
        np.array([p[2] for p in prims])) for l, prims in Bq[9]["shells"]]


def build_static(R):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    S = overlap.overlap([9, 9], xyz, charge=0, ao_order="oracle")
    cn = adapt.basis_cn([9, 9], xyz)[0]
    q_ = KB * math.sqrt(cn)
    kmap = {0: KS, 1: KP}
    sh = [{"l": l, "at": at, "exp": e * kmap[l], "coef": (c0 + c1 * q_).copy()}
          for at in (0, 1) for l, e, c0, c1 in _sh]
    Sh = overlap.overlap([9, 9], xyz, shells=sh)
    n = 8
    H0 = np.diag([-L2[LOFAO[i]] for i in range(n)]).astype(float)
    for i in range(n):
        for j in range(n):
            if ATOFAO[i] == ATOFAO[j]:
                continue
            li, lj = LOFAO[i], LOFAO[j]
            sig = (i % 4 in (0, 3)) and (j % 4 in (0, 3))
            if li == 0 and lj == 0:
                a = kW_s * kdiat_sg
            elif li != lj:
                a = (kW_s + kW_p) / 2 * kdiat_sg
            else:
                a = kW_p * (kdiat_sg if sig else kdiat_pi)
            H0[i, j] = a * (-(L2[li] + L2[lj]) / 2) * float(Sh[i, j])
    A = f2_stretch.acp_matrix([9, 9], xyz)
    obj = json.load(open(os.path.join(HERE, "data", "derived-constants.json")))[
        "hamiltonian_anchors"]["pzpz_object"]["closed_form_MEASURED"]
    OBJ = np.zeros((n, n))
    OBJ[3, 7] = OBJ[7, 3] = obj["c"] * math.erf(obj["a"] * R) / R
    # the EMPIRICAL short pieces (f2_shortfit; labeled working laws), with the
    # center-swap mirror rules (spz-type elements are antisymmetric under the swap)
    sf = os.path.join(HERE, "data", "f2-short.json")
    if INCLUDE_SHORT and os.path.exists(sf):
        fits = json.load(open(sf))["fits"]

        def ev(k, sval):
            f = fits[k]
            a_ = abs(sval)
            fns = {"s2,s4": (a_ * a_, a_ ** 4), "s,s2": (a_, a_ * a_),
                   "s2,s3": (a_ * a_, a_ ** 3), "s,s3": (a_, a_ ** 3)}[f["form"]]
            return f["c"][0] * fns[0] + f["c"][1] * fns[1]

        SH = np.zeros((n, n))
        for i, k in ((0, "d_s"), (4, "d_s"), (1, "d_px"), (2, "d_px"), (5, "d_px"),
                     (6, "d_px"), (3, "d_pz"), (7, "d_pz")):
            car = {"d_s": S[0, 4], "d_px": S[1, 5], "d_pz": S[3, 7]}[k]
            SH[i, i] = ev(k, car)
        SH[0, 4] = SH[4, 0] = ev("ss", S[0, 4])
        SH[3, 7] = SH[7, 3] = ev("pzpz", S[3, 7])
        for i, j in ((1, 5), (2, 6)):
            SH[i, j] = SH[j, i] = ev("pxpx", S[1, 5])
        v = ev("spz", S[0, 7])
        SH[0, 7] = SH[7, 0] = v
        SH[3, 4] = SH[4, 3] = -v
        v = ev("on_spz", S[0, 7])
        SH[0, 3] = SH[3, 0] = v
        SH[4, 7] = SH[7, 4] = -v
        OBJ = OBJ + SH
    return S, H0, A, OBJ


def gamma2_onsite():
    g = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            g[a, b] = S_F * 2 * U[a] * U[b] / (U[a] + U[b])
    return g


def fock(P, S, H0, A, OBJ, g2):
    n = 8
    Ps = P / 2.0
    m = np.diag(Ps @ S)
    # shell charges per atom (total, both spins)
    q = np.zeros((2, 2))                 # [atom][shell]
    for i in range(n):
        q[ATOFAO[i]][LOFAO[i]] += 2 * m[i]
    for at in range(2):
        q[at][0] = REF[0] - q[at][0]
        q[at][1] = REF[1] - q[at][1]
    v = np.zeros(n)
    for i in range(n):
        at, l = ATOFAO[i], LOFAO[i]
        v[i] += MU[l]                                          # ES1
        v[i] += g2[l][0] * q[at][0] + g2[l][1] * q[at][1]      # ES2 onsite
        v[i] += 2 * CX * L5[l] * m[i]                          # X same-shell per-AO
    F = H0 + A + OBJ - 0.5 * S * (v[:, None] + v[None, :])
    return F


def scf(R, iters=60):
    S, H0, A, OBJ = build_static(R)
    g2 = gamma2_onsite()
    n = 8
    s, Uo = np.linalg.eigh(S)
    X = Uo @ np.diag(1 / np.sqrt(s)) @ Uo.T
    P = np.zeros((n, n))
    for it in range(iters):
        F = fock(P, S, H0, A, OBJ, g2)
        w, Vp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Vp
        Pn = 2.0 * C[:, :7] @ C[:, :7].T
        if np.abs(Pn - P).max() < 1e-9:
            P = Pn
            break
        P = 0.6 * Pn + 0.4 * P
    F = fock(P, S, H0, A, OBJ, g2)
    w, Vp = np.linalg.eigh(X.T @ F @ X)
    return w, P


def main():
    for R in (3.8, 3.0, 2.668):
        w, P = scf(R)
        r = oracle.run([("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
        ref = sorted(x / K.EV for x in r["eps_ev"])
        d = [a - b for a, b in zip(sorted(w), ref)]
        worst = max(abs(x) for x in d)
        tag = "GATE" if R == 3.8 else "diag"
        print(f"R={R} [{tag}]  worst |d| = {worst:.4f} "
              f"{'PASS' if worst <= 0.02 and R == 3.8 else ''}")
        for a, b, dd in zip(sorted(w), ref, d):
            print(f"    {a:+.4f}  vs {b:+.4f}   d {dd:+.4f}")


if __name__ == "__main__":
    main()
