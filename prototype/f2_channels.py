"""f2_channels.py -- complete the F2 mu/CN channel measurements at ALL 10 stretch points.

The v2 forward gate's p-element failures peak at short R where the mu (L7) and CN (L3,
L1[8]) channel responses were INTERPOLATED from just 4 atlas points. This measures the five
responses at every stretch distance so the channels enter the forward model as data, not
interpolation. Writes data/f2-channels.json.
"""
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = [2.0, 2.3, 2.668, 3.0, 3.4, 3.8, 4.2, 4.6, 5.2, 6.0]
ELS = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}
D = 0.02
SLOTS = [("L7_s", (9, 7, 0)), ("L7_p", (9, 7, 1)), ("L3_s", (9, 3, 0)),
         ("L3_p", (9, 3, 1)), ("L1[8]k1cn", (9, 1, 8))]


def main():
    P = params.parse()
    eF = P["element"][9]
    vals = {"L7_s": eF["shells"][5][0], "L7_p": eF["shells"][5][1],
            "L3_s": eF["shells"][1][0], "L3_p": eF["shells"][1][1],
            "L1[8]k1cn": eF["l1"][8]}
    out = {}
    try:
        for R in RS:
            atoms = [("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)]
            rec = {}
            for name, key in SLOTS:
                params.write_perturbed({key: vals[name] + D})
                Fp = fock_recon.fock_ao(atoms)["F"]
                params.write_perturbed({key: vals[name] - D})
                Fm = fock_recon.fock_ao(atoms)["F"]
                rec[name] = {el: float((Fp[ix] - Fm[ix]) / (2 * D))
                             for el, ix in ELS.items()}
            out[str(R)] = rec
            print(f"  R={R:5.2f} done (L7_p on pxpx {rec['L7_p']['pxpx']:+.5f})", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "f2-channels.json"), "w"), indent=1)
    print("wrote data/f2-channels.json")


if __name__ == "__main__":
    main()
