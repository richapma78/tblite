"""h2plus_fock.py -- measure the exchange Fock DIRECTLY on the 1-electron manifold.

H2+ (UKS, one alpha electron): the beta channel is EMPTY, and for a P-quadratic exchange
energy the beta-Fock's exchange is exactly zero (linear gradient at P^beta = 0). Both spin
channels share H0, ES, ACP (charge-driven, spin-blind), so

    F^alpha - F^beta = X^alpha(P^alpha) + [spin-polarization term]

and the spin term is MEASURED by its own slot (L1[9] = W): spin-content = W * dF/dW per
element. Everything else cancels -- no forward-model stacking at all. The four printed
eigenvalues (alpha pair + beta pair) invert with our S (charge-adapted basis at q = +1).

DISCRIMINATOR (pre-declared): the Mulliken-POPULATION form predicts X(H2+) = X(H2-per-spin)
(both have m = 1/2); any P-MATRIX-linear form predicts X(H2+) = X(H2)/2 (P^alpha is half of
H2's per-spin P). The measured curve picks. Writes data/h2plus-fock.json.
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
    r = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)], charge=1, uhf=1)
    eps = [x / K.EV for x in r["eps_ev"]]
    assert len(eps) == 4, eps
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    S12 = float(overlap.overlap([1, 1], xyz, charge=1)[0, 1])
    Fa = invert(eps[0:2], S12)
    Fb = invert(eps[2:4], S12)
    return Fa, Fb, S12, r


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
            out.append({"R": R, "S12": S12, "Fa": Fa, "Fb": Fb,
                        "spin11": spin11, "spin12": spin12,
                        "X11": X11, "X12": X12, "Ex": r["terms"]["Ex (Mulliken)"],
                        "q_H": r["eeq"][0]["q"]})
            print(f"  R={R:4.2f}  X11 {X11:+.5f}  X12 {X12:+.5f}  "
                  f"(spin {spin11:+.4f}/{spin12:+.4f})  q_H {r['eeq'][0]['q']:+.3f}",
                  flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "h2plus-fock.json"), "w"), indent=1)
    # discriminator table vs the H2 X_required (x-placement.json)
    xp = {p["R"]: p for p in json.load(open(os.path.join(HERE, "data",
                                                         "x-placement.json")))["points"]}
    print("\nDISCRIMINATOR: X11(H2+) vs X11_req(H2) -- equal => population form; "
          "half => P-matrix form:")
    for rec in out:
        R = rec["R"]
        if R in xp:
            print(f"  R={R:4.2f}  H2+ {rec['X11']:+.5f}   H2 {xp[R]['X11_req']:+.5f}   "
                  f"ratio {rec['X11'] / xp[R]['X11_req']:+.3f}")


if __name__ == "__main__":
    main()
