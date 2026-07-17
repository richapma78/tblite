"""h2triplet_fock.py -- the FOURTH exchange manifold: triplet H2 (.UHF 2, beta empty).

Same instrument as h2plus_fock (F^alpha - F^beta = X + spin; spin via the W-slot FD), but a
radically different alpha density: BOTH MOs occupied, so in closed form

    m^alpha_A = 1 (vs 1/2 in H2+ and per-spin singlet)
    P^alpha_12 = -S/(1-S^2)   (NEGATIVE and S-shaped, vs +1/(2(1+S)))

The diagonal tests the m-scaling of -gamma_on*m directly; the off-diagonal's sign flip
separates the S-carried from the P12-carried pieces that were degenerate on the singlet-type
manifolds. Writes data/h2triplet-fock.json.
"""
import json
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

PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = [1.0, 1.4, 2.0, 2.5, 3.0, 4.0, 5.0]
DW = 0.02


def invert(epair, S12):
    eb, ea = sorted(epair)
    F11 = ((1 + S12) * eb + (1 - S12) * ea) / 2
    F12 = ((1 + S12) * eb - (1 - S12) * ea) / 2
    return F11, F12


def fock_pair(R):
    r = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)], charge=0, uhf=2)
    eps = [x / K.EV for x in r["eps_ev"]]
    assert len(eps) == 4, eps
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    S12 = float(overlap.overlap([1, 1], xyz, charge=0)[0, 1])
    return invert(eps[0:2], S12), invert(eps[2:4], S12), S12, r


def main():
    e = params.parse()["element"][1]
    W = e["l1"][9]
    out = []
    try:
        for R in RS:
            params.write_perturbed({})
            Fa, Fb, S12, r = fock_pair(R)
            params.write_perturbed({(1, 1, 9): W + DW})
            Fap, Fbp, _s, _r = fock_pair(R)
            params.write_perturbed({(1, 1, 9): W - DW})
            Fam, Fbm, _s, _r = fock_pair(R)
            spin11 = W * ((Fap[0] - Fbp[0]) - (Fam[0] - Fbm[0])) / (2 * DW)
            spin12 = W * ((Fap[1] - Fbp[1]) - (Fam[1] - Fbm[1])) / (2 * DW)
            X11 = (Fa[0] - Fb[0]) - spin11
            X12 = (Fa[1] - Fb[1]) - spin12
            out.append({"R": R, "S12": S12, "X11": X11, "X12": X12,
                        "spin11": spin11, "spin12": spin12,
                        "Ex": r["terms"].get("Ex (Mulliken)"),
                        "gap_ev": r["gap_ev"]})
            print(f"  R={R:4.2f}  X11 {X11:+.5f}  X12 {X12:+.5f}  "
                  f"(spin {spin11:+.4f}/{spin12:+.4f})  Ex {r['terms'].get('Ex (Mulliken)'):+.5f}",
                  flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "h2triplet-fock.json"), "w"), indent=1)
    print("wrote data/h2triplet-fock.json")
    print("\npredictions to test offline: X11 ~ -gamma_on*1 (m=1); X12's P12-part flips sign")


if __name__ == "__main__":
    main()
