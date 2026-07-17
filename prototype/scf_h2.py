"""scf_h2.py -- the FIRST SCF ASSEMBLY: H2 from gated parts, no oracle in the loop.

Pieces (every one previously gated or measured):
  S        our overlap (machine-exact vs the binary)
  H0       the forward-gated assembly: off-diag amp*(-L2)*Pi(b)*S~(k, kb) (worst 0.0094 on
           13 points); diag -L2 (the L3*CN term bounded <= 0.005 at R <= 1.4)
  ES1      F -= 1/2 * (mu_A corr_A + mu_B corr_B) * S_munu -- the decoded -mu*S law with the
           MEASURED corr(R) = -C_mu(R) (data/mu-cmu-h2.json; exactly 1.000 at R >= 4)
  ES2      gamma2 = s(H)*U (onsite, decoded) + offsite scan-level kernel; enters via
           F += -1/2 S (v_A + v_B), v = gamma*q -- ZERO at the symmetric fixed point, kept
           for iteration stability only
  X        Eq. 153 with the MEASURED kernel: gamma_on = 0.360130 (atom-pinned),
           gamma_off(R) interpolated from the 17-point curve (data/ex-kernel + ex-tail)
  ACP      the analytic projector matrix (validated against FD both on- and off-diagonal)

GATES (pre-declared, at R = 1.4, 2.5, 4.0):
  G-eps    both eigenvalues within 0.02 Eh of the printed pair
  G-Ex     Ex = 1/2 Tr(F^X P) within 0.002 Eh of the printed 'Ex (Mulliken)'
  G-ES1    E1 = sum mu*corr*q ~ 0 at q = 0: |printed ES1 - ours| <= 0.002
  (the full electronic-energy gate needs the remaining energy bookkeeping -- next push)
Failures are reported per piece; nothing is tuned here.
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
eH = P_["element"][1]
AMP = P_["globals"]["g1"][0] * eH["l8"][0]
L2, MU, L5 = eH["shells"][0][0], eH["shells"][5][0], eH["shells"][3][0]
K_FIT, KB_FIT, B_SHP = 1.110, 0.179, 0.0315
GAM_ON = 0.360130

_ker = json.load(open(os.path.join(HERE, "data", "ex-kernel.json")))
_tail = json.load(open(os.path.join(HERE, "data", "ex-tail.json")))
_kR = [r["R"] for r in _ker["rows"]] + [r["R"] for r in _tail["tail"]]
_kV = [r["gamma_off"] for r in _ker["rows"]] + [r["gamma_off"] for r in _tail["tail"]]
_cmu = json.load(open(os.path.join(HERE, "data", "mu-cmu-h2.json")))
_cR = sorted(float(x) for x in _cmu)
_cV = [_cmu[str(x)]["C_mu"] for x in _cR]

Bq = basisq.parse()
_prims = Bq[1]["shells"][0][1]
_e0 = np.array([p[0] for p in _prims])
_c0 = np.array([p[1] for p in _prims])
_c1 = np.array([p[2] for p in _prims])


def build(R):
    """All P-independent matrices for H2 at R (Bohr)."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    S = overlap.overlap([1, 1], xyz)
    cn = adapt.basis_cn([1, 1], xyz)[0]
    q_ = KB_FIT * math.sqrt(cn)
    sh = [{"l": 0, "at": i, "exp": _e0 * K_FIT, "coef": (_c0 + _c1 * q_).copy()}
          for i in (0, 1)]
    Sham = overlap.overlap([1, 1], xyz, shells=sh)
    H0 = np.diag([-L2, -L2]).astype(float)
    H0[0, 1] = H0[1, 0] = AMP * (-L2) * (1 + B_SHP * R) ** 2 * float(Sham[0, 1])
    A = f2_stretch.acp_matrix([1, 1], xyz)
    corr = -float(np.interp(R, _cR, _cV))
    ES1 = -0.5 * (MU * corr + MU * corr) * S          # both centers H
    gam = np.array([[GAM_ON, float(np.interp(R, _kR, _kV))],
                    [float(np.interp(R, _kR, _kV)), GAM_ON]])
    return S, H0, A, ES1, gam


def ex_energy(P, S, gam):
    """The GATED energy form (ex_extract conventions): spin-resolved Eq. 151 with
    P^sigma = P/2, doubled -- its VALUES matched the printed Ex on every tested
    configuration (L5 grids, R grids, atom)."""
    Ph = P / 2.0
    n = P.shape[0]
    E = 0.0
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                for kap in range(n):
                    E -= Ph[mu, nu] * S[mu, lam] * Ph[lam, kap] * S[kap, nu] * (
                        gam[mu, kap] + gam[mu, nu] + gam[lam, kap] + gam[lam, nu]) / 16.0
    return 2.0 * E


def fock_x(P, S, gam):
    """The MEASURED exchange-Fock skeleton (forty-seventh push): the onsite-kernel part is
    IDENTIFIED as the Mulliken potential of v_A = gamma_on * m_A^sigma (own same-spin
    population): F = -1/2 S o (v_mu + v_nu). This unifies the exact atom anchor, the
    diagonal m-scaling ladder, and both molecules' off-diagonal asymptotics. The smaller
    gamma_off-carried remainder (singlet long-R: -2*c_x*gamma_off*P12, flagged) is NOT yet
    identified and NOT included -- the gate reports what the identified part alone gives.
    NOTE: deliberately non-variational, as the binary's own Fock is."""
    Ps = P / 2.0
    m = np.diag(Ps @ S)                                # per-spin Mulliken populations
    v = np.diag(gam) * m                               # onsite kernel x own population
    return -0.5 * S * (v[:, None] + v[None, :])


def scf(R, iters=40):
    S, H0, A, ES1, gam = build(R)
    n = 2
    P = np.zeros((n, n))
    e, C = None, None
    for _ in range(iters):
        FX = fock_x(P, S, gam) if P.any() else np.zeros((n, n))
        F = H0 + ES1 + A + FX
        # generalized eigenproblem via Loewdin
        s, U = np.linalg.eigh(S)
        X = U @ np.diag(1 / np.sqrt(s)) @ U.T
        w, Vp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Vp
        Pn = 2.0 * np.outer(C[:, 0], C[:, 0])
        if np.abs(Pn - P).max() < 1e-10:
            P = Pn
            break
        P = 0.5 * (Pn + P)
    FX = fock_x(P, S, gam)
    F = H0 + ES1 + A + FX
    s, U = np.linalg.eigh(S)
    X = U @ np.diag(1 / np.sqrt(s)) @ U.T
    w, Vp = np.linalg.eigh(X.T @ F @ X)
    Ex = ex_energy(P, S, gam)
    return {"eps": w, "P": P, "Ex": Ex, "F": F, "FX": FX}


def main():
    print("H2 SCF from gated parts vs the oracle:")
    ok = True
    for R in (1.4, 2.5, 4.0):
        res = scf(R)
        r = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
        eps_ref = sorted(x / K.EV for x in r["eps_ev"])
        d0 = res["eps"][0] - eps_ref[0]
        d1 = res["eps"][1] - eps_ref[1]
        dEx = res["Ex"] - r["terms"]["Ex (Mulliken)"]
        dE1 = 0.0 - r["terms"]["ES1 (charge SIE)"]
        g_eps = abs(d0) <= 0.02 and abs(d1) <= 0.02
        g_ex = abs(dEx) <= 0.002
        g_e1 = abs(dE1) <= 0.002
        ok &= g_eps and g_ex and g_e1
        print(f"  R={R:4.2f}  eps ({res['eps'][0]:+.4f}, {res['eps'][1]:+.4f}) vs "
              f"({eps_ref[0]:+.4f}, {eps_ref[1]:+.4f})  d = ({d0:+.4f}, {d1:+.4f}) "
              f"[{'OK' if g_eps else 'MISS'}]")
        print(f"          Ex {res['Ex']:+.5f} vs {r['terms']['Ex (Mulliken)']:+.5f} "
              f"(d {dEx:+.5f}) [{'OK' if g_ex else 'MISS'}]   ES1 ref "
              f"{r['terms']['ES1 (charge SIE)']:+.5f} (ours 0) [{'OK' if g_e1 else 'MISS'}]",
              flush=True)
    print("SCF ASSEMBLY GATE " + ("PASSED" if ok else "FAILED (report per piece above)"))


if __name__ == "__main__":
    main()
