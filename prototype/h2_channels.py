"""h2_channels.py -- complete the H2 mu/L3 channel measurements: both Fock elements, all
13 stretch points.

Why: the ES1 term is mu-LINEAR (F contains -mu*(1 + k1cn*sqrt(CN_mu))*metric), so the
mu-slot response ALONE carries the complete channel (the k1cn structure rides inside mu's
coefficient) -- mu * dF/dmu IS the term, on the diagonal as well as off. With both elements
measured everywhere, the exchange extraction gets a clean path:

    X11_required = F11(reconstructed) - (-L2) - mu*resp11_mu - L3*resp11_L3 - ACP11

which involves NO forward-model fits at all (H0's diagonal is exactly -L2 + the measured
L3 content) -- the 1e-3-grade extraction the placement fit needs. Writes
data/h2-channels.json.
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
RS = [0.9, 1.0, 1.2, 1.4, 1.7, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
D = 0.02
SLOTS = [("L7", (1, 7, 0)), ("L3", (1, 3, 0)), ("L2", (1, 2, 0))]


def main():
    e = params.parse()["element"][1]
    vals = {"L7": e["shells"][5][0], "L3": e["shells"][1][0], "L2": e["shells"][0][0]}
    out = {}
    try:
        for R in RS:
            atoms = [("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)]
            rec = {}
            for name, key in SLOTS:
                params.write_perturbed({key: vals[name] + D})
                Fp = fock_recon.fock_ao(atoms)["F"]
                params.write_perturbed({key: vals[name] - D})
                Fm = fock_recon.fock_ao(atoms)["F"]
                rec[name] = {"d11": float((Fp[0, 0] - Fm[0, 0]) / (2 * D)),
                             "d12": float((Fp[0, 1] - Fm[0, 1]) / (2 * D))}
            out[str(R)] = rec
            print(f"  R={R:4.2f}  mu(d11 {rec['L7']['d11']:+.4f}, d12 {rec['L7']['d12']:+.4f})"
                  f"  L3(d11 {rec['L3']['d11']:+.4f}, d12 {rec['L3']['d12']:+.4f})"
                  f"  L2(d11 {rec['L2']['d11']:+.4f})", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "h2-channels.json"), "w"), indent=1)
    print("wrote data/h2-channels.json")


if __name__ == "__main__":
    main()
