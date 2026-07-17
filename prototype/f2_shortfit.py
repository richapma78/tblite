"""f2_shortfit.py -- extract and fit the F2 short pieces (the H2 route, p-shell edition).

Per stretch point: remainder = F_recon - [H0 + ACP + ES1 + ES2 + X-skeleton + object]
evaluated at the oracle's density. Fit each element class's curve to 2-parameter forms in
the overlap/bond objects; install as LABELED EMPIRICAL closed forms; re-gate f2_scf at r_e.
(The k_p metric puzzle's misfit is absorbed here and stays flagged.)
Writes data/f2-short.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import f2_scf  # noqa: E402
f2_scf.INCLUDE_SHORT = False
import fock_recon  # noqa: E402

RS = [2.0, 2.3, 2.668, 3.0, 3.4, 3.8, 4.2, 4.6, 5.2, 6.0]
CLASSES = {"d_s": (0, 0), "d_px": (1, 1), "d_pz": (3, 3),
           "ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5),
           "on_spz": (0, 3)}


def main():
    g2 = f2_scf.gamma2_onsite()
    data = {k: [] for k in CLASSES}
    svals = []
    for R in RS:
        rec = fock_recon.fock_ao([("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
        S_, H0, A, OBJ = f2_scf.build_static(R)
        F_asm = f2_scf.fock(rec["state"]["P"], S_, H0, A, OBJ, g2)
        rem = rec["F"] - F_asm
        svals.append({"R": R, "s_ss": float(S_[0, 4]), "s_spz": float(S_[0, 7]),
                      "s_pzpz": float(S_[3, 7]), "s_pxpx": float(S_[1, 5])})
        for k, ix in CLASSES.items():
            data[k].append(float(rem[ix]))
        print(f"  R={R:5.2f}  " + "  ".join(f"{k} {rem[ix]:+.4f}" for k, ix in
                                            list(CLASSES.items())[:5]), flush=True)
    Rv = np.array(RS)

    # fit each class: c1*f1 + c2*f2 with the relevant overlap as the carrier
    fits = {}
    print("\nfits (2-param, rms):")
    for k in CLASSES:
        y = np.array(data[k])
        skey = {"d_s": "s_ss", "d_px": "s_pxpx", "d_pz": "s_pzpz", "ss": "s_ss",
                "spz": "s_spz", "pzpz": "s_pzpz", "pxpx": "s_pxpx",
                "on_spz": "s_spz"}[k]
        sv = np.array([abs(r[skey]) for r in svals])
        best = None
        for name, f1, f2 in (("s2,s4", sv * sv, sv ** 4),
                             ("s,s2", sv, sv * sv),
                             ("s2,s3", sv * sv, sv ** 3),
                             ("s,s3", sv, sv ** 3)):
            Am = np.vstack([f1, f2]).T
            c, *_ = np.linalg.lstsq(Am, y, rcond=None)
            rms = float(np.sqrt(np.mean((Am @ c - y) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, name, c)
        rms, name, c = best
        fits[k] = {"form": name, "c": [float(x) for x in c], "rms": rms,
                   "carrier": skey}
        print(f"  {k:6s} [{name}] c = ({c[0]:+.4f}, {c[1]:+.4f})  rms {rms:.2e}")
    json.dump({"remainders": data, "svals": svals, "fits": fits},
              open(os.path.join(HERE, "data", "f2-short.json"), "w"), indent=1)
    print("wrote data/f2-short.json")


if __name__ == "__main__":
    main()
