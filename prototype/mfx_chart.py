"""mfx_chart.py -- chart the binary's MFX exchange kernel gamma(R; U) empirically.

Eq. 149 (SI) puts the shell Hubbards U^MFX (= L5, the only shell-wise row left) in the
KERNEL DENOMINATOR via favg -- onsite gamma(0) = alpha/U, so |Ex| would SHRINK with L5.
The gated atom law says the opposite (E_x = -c_x*L5*n^2, |Ex| GROWING linearly). One of
them is wrong about the binary: SET L5(H) on a grid and read the printed Ex (atom, and H2
at three R) plus the reconstructed F12/F11. The chart's shape is the verdict:
    linear-increasing |Ex(L5)|  -> the binary's kernel carries U in the NUMERATOR
                                   (paper-vs-binary divergence D4, to be registered)
    ~1/L5                       -> Eq. 149 as printed; the atom law needs re-reading
Also charts X12(L5) offsite at fixed R (the R-dependence separates the erf/screening
factors for the later (alpha, omega, k1, k2) fit). Writes data/mfx-chart.json.
"""
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
GRID = [0.5, 1.0, 2.0, 3.6548180166, 5.0, 8.0]


def main():
    out = {"atom": [], "h2": {}}
    try:
        for v in GRID:
            params.write_perturbed({(1, 5, 0): v})
            r = oracle.run([("H", 0, 0, 0)], uhf=1)
            ex = r["terms"].get("Ex (Mulliken)")
            out["atom"].append({"L5": v, "Ex": ex})
            print(f"  atom L5={v:6.3f}  Ex {ex:+.6f}", flush=True)
        for R in (1.4, 2.5, 4.0):
            rows = []
            for v in GRID:
                params.write_perturbed({(1, 5, 0): v})
                rr = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
                ex = rr["terms"].get("Ex (Mulliken)")
                fr = fock_recon.fock_ao([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
                rows.append({"L5": v, "Ex": ex, "F12": float(fr["F"][0, 1]),
                             "F11": float(fr["F"][0, 0])})
                print(f"  H2 R={R} L5={v:6.3f}  Ex {ex:+.6f}  F12 {fr['F'][0, 1]:+.6f}",
                      flush=True)
            out["h2"][str(R)] = rows
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump(out, open(os.path.join(HERE, "data", "mfx-chart.json"), "w"), indent=1)

    print("\nVERDICT SHAPE (atom): Ex vs L5 -- linear or 1/L5?")
    a = out["atom"]
    for i in range(1, len(a)):
        dE = a[i]["Ex"] - a[i - 1]["Ex"]
        dL = a[i]["L5"] - a[i - 1]["L5"]
        print(f"  L5 {a[i - 1]['L5']:.2f}->{a[i]['L5']:.2f}: dEx/dL5 = {dE / dL:+.6f}")


if __name__ == "__main__":
    main()
