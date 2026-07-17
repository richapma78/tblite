"""ex_extract.py -- extract the binary's exchange kernel gamma(R) from printed Ex + restart P.

D4 (mfx_chart, banked): the binary's Ex is EXACTLY linear in L5 (dEx/dL5 = -c_x to 6 digits,
c_x = s/9.59) -- U sits in the kernel NUMERATOR, refuting Eq. 149's denominator placement.
With that, Eq. 151's Mulliken structure for H2 (one s AO per atom) has exactly TWO kernel
values: gamma_on (R = 0, same atom) and gamma_off(R). Both P and S are known per point
(restart + our overlap), so the printed Ex INVERTS:

    Ex(R) = A(P, S) * gamma_on + B(P, S) * gamma_off(R)

with A, B computed from the exact Eq-151 contraction (both spins; RKS closed shell). The
atom (single alpha electron) pins gamma_on independently; the stretch then yields the
gamma_off(R) CURVE -- the object the (alpha, omega, k1, k2) fit must reproduce.
Writes data/ex-kernel.json.
"""
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = [0.9, 1.0, 1.2, 1.4, 1.7, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]


def ex_terms(P, S, gam):
    """Eq. 151 spin-summed for a CLOSED SHELL with total-P convention: E = -1/16 * sum
    P_munu S_mulam P_lamkap S_kapnu (g_mukap + g_munu + g_lamkap + g_lamnu).
    (RKS: P_sigma = P/2 per spin, two spins -> 2 * (1/8) * (1/4) = 1/16.) The overall
    prefactor convention is FITTED below anyway (A, B carry it); this fixes shape only."""
    n = P.shape[0]
    E = 0.0
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                for kap in range(n):
                    E -= P[mu, nu] * S[mu, lam] * P[lam, kap] * S[kap, nu] * (
                        gam[mu, kap] + gam[mu, nu] + gam[lam, kap] + gam[lam, nu]) / 16.0
    return E


def main():
    out = {"rows": []}
    try:
        params.write_perturbed({})
        # ---- the atom pins gamma_on: 1 alpha electron, P_aa = 1 (per spin)
        r = oracle.run([("H", 0, 0, 0)], uhf=1)
        ex_atom = r["terms"]["Ex (Mulliken)"]
        # spin-resolved: E = -1/8 * sum P^s P^s S S (4 gammas); atom: P^alpha = 1, S = 1
        # -> E = -1/8 * 1*1*1*1*4*gamma_on = -gamma_on/2
        gam_on = -2.0 * ex_atom
        print(f"atom: Ex {ex_atom:+.6f} -> gamma_on = {gam_on:+.6f}  "
              f"(c_x*L5 = {0.049268 * 3.6548180166:+.6f})", flush=True)
        out["gamma_on"] = gam_on
        for R in RS:
            atoms = [("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)]
            rr = oracle.run(atoms)
            ex = rr["terms"]["Ex (Mulliken)"]
            st = fock_recon.fock_ao(atoms)
            P, S = st["state"]["P"], st["S"]
            # closed-shell spin-resolved: P^sigma = P/2, E = 2 * E_sigma(P/2)
            Ph = P / 2.0
            g_on = np.array([[1.0, 0.0], [0.0, 1.0]])
            g_off = np.array([[0.0, 1.0], [1.0, 0.0]])
            A = 2 * ex_terms(Ph, S, g_on)
            B = 2 * ex_terms(Ph, S, g_off)
            gam_off = (ex - A * gam_on) / B
            out["rows"].append({"R": R, "Ex": ex, "P12": float(P[0, 1]), "S12": float(S[0, 1]),
                                "A": A, "B": B, "gamma_off": float(gam_off)})
            print(f"  R={R:4.2f}  Ex {ex:+.6f}  P12 {P[0, 1]:+.4f}  A {A:+.4f}  B {B:+.4f}"
                  f"  gamma_off {gam_off:+.6f}", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "ex-kernel.json"), "w"), indent=1)
    print("\ngamma_off(R) shape diagnostics:")
    rows = out["rows"]
    for a, b in zip(rows, rows[1:]):
        if abs(a["gamma_off"]) > 1e-8 and abs(b["gamma_off"]) > 1e-8:
            print(f"  {a['R']:4.2f}->{b['R']:4.2f}: ratio {a['gamma_off'] / b['gamma_off']:.3f}"
                  f"   1/R ratio {b['R'] / a['R']:.3f}")


if __name__ == "__main__":
    main()
