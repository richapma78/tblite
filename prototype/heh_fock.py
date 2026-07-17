"""heh_fock.py -- the HETERONUCLEAR exchange decomposition: HeH2+ (one electron, beta empty).

Two elements' kernels in one clean spin-subtracted system. Per R:
  - four printed eigenvalues -> alpha/beta 2x2 inversions with our S (charge +2 basis chain)
  - spin term measured via BOTH W slots (He and H)
  - P^alpha reconstructed in closed form from the PRINTED Mulliken populations
    (m_He + m_H = 1; c1*c2 solved from the overlap-population quadratic)
  - the cross kernel gamma_HeH(R) extracted from the printed Ex via the Eq-151 contraction
    (basis gamma-matrices: He-only, H-only, cross-only -> A, B, C -> one unknown per point)
  - LAW TESTS with unequal populations: diag X_AA vs -gamma_on(A)*m_A (the population law);
    off-diag vs -1/2 s (v_He + v_H); the (PS - I) branch.
gamma_on(He) = 0.899150 (calibrated, slope-exact), gamma_on(H) = 0.360130.
Writes data/heh-fock.json.
"""
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402
from scf_h2 import ex_energy  # noqa: E402

PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = [1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0]
GHE, GH = 0.899150, 0.360130
DW = 0.02


def invert(epair, s):
    eb, ea = sorted(epair)
    return (((1 + s) * eb + (1 - s) * ea) / 2, ((1 + s) * eb - (1 - s) * ea) / 2)


def run_pair(R):
    r = oracle.run([("He", 0, 0, 0), ("H", 0, 0, R / K.BOHR)], charge=2, uhf=1)
    eps = [x / K.EV for x in r["eps_ev"]]
    assert len(eps) == 4, eps
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    s = float(overlap.overlap([2, 1], xyz, charge=2)[0, 1])
    return invert(eps[0:2], s), invert(eps[2:4], s), s, r


def p_alpha(mHe, mH, s):
    """1-electron 2-AO density from Mulliken populations (bonding sign)."""
    inv = 1.0 / (s * s) - 1.0
    t = (-1.0 + math.sqrt(1.0 + 4.0 * inv * mHe * mH)) / (2.0 * inv)
    c1 = math.sqrt(max(mHe - t, 0.0))
    c2 = math.sqrt(max(mH - t, 0.0))
    return np.array([[c1 * c1, c1 * c2], [c1 * c2, c2 * c2]])


def main():
    P = params.parse()
    WHe, WH = P["element"][2]["l1"][9], P["element"][1]["l1"][9]
    out = []
    try:
        for R in RS:
            params.write_perturbed({})
            Fa, Fb, s, r = run_pair(R)
            spins = {}
            for tag, key, w0 in (("He", (2, 1, 9), WHe), ("H", (1, 1, 9), WH)):
                params.write_perturbed({key: w0 + DW})
                Fap, Fbp, _s, _r = run_pair(R)
                params.write_perturbed({key: w0 - DW})
                Fam, Fbm, _s, _r = run_pair(R)
                d11 = w0 * ((Fap[0] - Fbp[0]) - (Fam[0] - Fbm[0])) / (2 * DW)
                d12 = w0 * ((Fap[1] - Fbp[1]) - (Fam[1] - Fbm[1])) / (2 * DW)
                spins[tag] = (d11, d12)
            spin11 = spins["He"][0] + spins["H"][0]
            spin12 = spins["He"][1] + spins["H"][1]
            X11 = (Fa[0] - Fb[0]) - spin11
            X12 = (Fa[1] - Fb[1]) - spin12
            mHe, mH = r["pops"][0]["s"], r["pops"][1]["s"]
            Pa = p_alpha(mHe, mH, s)
            S = np.array([[1.0, s], [s, 1.0]])
            # kernel decomposition of the printed Ex (per-spin form: pass 2*P^alpha)
            gHe_m = np.array([[1.0, 0.0], [0.0, 0.0]])
            gH_m = np.array([[0.0, 0.0], [0.0, 1.0]])
            gx_m = np.array([[0.0, 1.0], [1.0, 0.0]])
            A_ = ex_energy(2 * Pa, S, gHe_m)
            B_ = ex_energy(2 * Pa, S, gH_m)
            C_ = ex_energy(2 * Pa, S, gx_m)
            ex = r["terms"]["Ex (Mulliken)"]
            gx = (ex - A_ * GHE - B_ * GH) / C_
            rec = {"R": R, "s": s, "X11": X11, "X12": X12, "mHe": mHe, "mH": mH,
                   "P12": float(Pa[0, 1]), "Ex": ex, "gamma_cross": float(gx),
                   "spin11": spin11, "spin12": spin12}
            out.append(rec)
            print(f"  R={R:4.2f}  X11 {X11:+.5f}  X12 {X12:+.5f}  mHe {mHe:.4f}  "
                  f"g_cross {gx:+.5f}", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "heh-fock.json"), "w"), indent=1)

    print("\nLAW TESTS (heteronuclear, unequal populations):")
    print("   R     X11+gHe*mHe   X12+s/2*(vHe+vH)   (PS-I)_12   g_cross")
    for rec in out:
        v = GHE * rec["mHe"] + GH * rec["mH"]
        d11 = rec["X11"] + GHE * rec["mHe"]
        d12 = rec["X12"] + 0.5 * rec["s"] * v
        Pa = p_alpha(rec["mHe"], rec["mH"], rec["s"])
        S = np.array([[1.0, rec["s"]], [rec["s"], 1.0]])
        psi = (Pa @ S - np.eye(2))[0, 1]
        print(f"  {rec['R']:4.2f}   {d11:+.5f}      {d12:+.5f}        {psi:+.4f}    "
              f"{rec['gamma_cross']:+.4f}")


if __name__ == "__main__":
    main()
