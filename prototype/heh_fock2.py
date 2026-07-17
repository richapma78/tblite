"""heh_fock2.py -- the heteronuclear exchange decomposition, PROPERLY (UKS reconstruction).

v1's symmetric inversion was invalid heteronuclear (bug caught + flagged). v2 uses the
per-spin full reconstruction F_sigma = S C_sigma eps_sigma C_sigma^T S from the decoded UKS
restart -- no symmetry assumption.

GATE first (pre-declared): H2+ via the UKS path must reproduce the old symmetric-inversion
X11/X12 (valid by symmetry) within 1e-4 at R = 1.4 and 2.5. Then HeH2+:
    X_el = (F_a - F_b)_el - spin-FD_el   (both W slots)
LAW TESTS with unequal populations (m_He ~ 0.98, m_H ~ 0.02):
    diag He: X_HeHe vs -gamma_on(He)*m_He      diag H: X_HH vs -gamma_on(H)*m_H
    off-diag vs -1/2 s (v_He + v_H) and the (PS - I) branch.
Writes data/heh-fock2.json.
"""
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402

PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = [1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0]
GHE, GH = 0.899150, 0.360130
DW = 0.02


def dF(atoms, charge):
    r = fock_recon.fock_ao(atoms, charge=charge, uhf=1)
    return r["F_a"] - r["F_b"], r


def main():
    P = params.parse()
    WHe, WH = P["element"][2]["l1"][9], P["element"][1]["l1"][9]
    try:
        # ---- GATE: H2+ UKS path vs the banked symmetric-inversion values
        params.write_perturbed({})
        hp = json.load(open(os.path.join(HERE, "data", "h2plus-fock.json")))
        ok = True
        for rec in hp:
            if rec["R"] not in (1.4, 2.5):
                continue
            atoms = [("H", 0, 0, 0), ("H", 0, 0, rec["R"] / K.BOHR)]
            D, r = dF(atoms, 1)
            # the banked X = dF - spin; compare the RAW dF against (X + spin) banked
            ref11 = rec["X11"] + rec["spin11"]
            ref12 = rec["X12"] + rec["spin12"]
            d11 = abs(D[0, 0] - ref11)
            d12 = abs(D[0, 1] - ref12)
            ok &= d11 < 1e-4 and d12 < 1e-4
            print(f"  GATE H2+ R={rec['R']}: dF11 {D[0, 0]:+.6f} (ref {ref11:+.6f}, "
                  f"d {d11:.1e})  dF12 {D[0, 1]:+.6f} (ref {ref12:+.6f}, d {d12:.1e})  "
                  f"ortho {r['ortho_err']:.1e}", flush=True)
        print("  UKS-vs-symmetric GATE " + ("PASSED" if ok else "FAILED -- stop"), flush=True)
        if not ok:
            return

        # ---- HeH2+ proper
        out = []
        for R in RS:
            params.write_perturbed({})
            atoms = [("He", 0, 0, 0), ("H", 0, 0, R / K.BOHR)]
            D0, r0 = dF(atoms, 2)
            spins = np.zeros((2, 2))
            for key, w0 in (((2, 1, 9), WHe), ((1, 1, 9), WH)):
                params.write_perturbed({key: w0 + DW})
                Dp, _ = dF(atoms, 2)
                params.write_perturbed({key: w0 - DW})
                Dm, _ = dF(atoms, 2)
                spins += w0 * (Dp - Dm) / (2 * DW)
            X = D0 - spins
            Pa = r0["state"]["P_a"]
            S = r0["S"]
            m = np.diag(Pa @ S)
            rec = {"R": R, "s": float(S[0, 1]),
                   "X_HeHe": float(X[0, 0]), "X_HH": float(X[1, 1]), "X_HeH": float(X[0, 1]),
                   "m_He": float(m[0]), "m_H": float(m[1]), "P12": float(Pa[0, 1]),
                   "PS_I_12": float((Pa @ S - np.eye(2))[0, 1]),
                   "ortho": r0["ortho_err"], "dens": r0["dens_err"]}
            out.append(rec)
            print(f"  R={R:4.2f}  X_HeHe {X[0, 0]:+.5f}  X_HH {X[1, 1]:+.5f}  "
                  f"X_HeH {X[0, 1]:+.5f}  m ({m[0]:.4f}, {m[1]:.4f})", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "heh-fock2.json"), "w"), indent=1)

    print("\nLAW TESTS:")
    print("   R    X_HeHe+gHe*mHe   X_HH+gH*mH   X_HeH+s/2*(vHe+vH)   (PS-I)12")
    for rec in out:
        v = GHE * rec["m_He"] + GH * rec["m_H"]
        print(f"  {rec['R']:4.2f}   {rec['X_HeHe'] + GHE * rec['m_He']:+.5f}       "
              f"{rec['X_HH'] + GH * rec['m_H']:+.5f}     "
              f"{rec['X_HeH'] + 0.5 * rec['s'] * v:+.5f}         {rec['PS_I_12']:+.4f}")


if __name__ == "__main__":
    main()
