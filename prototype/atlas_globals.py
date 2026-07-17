"""atlas_globals.py -- offdiag atlas, part 2: the GLOBAL rows + the L8 quartet, at short R.

The SI (Sec 1.7) says the EHT carries 8 globals -- k^W_l (Wolfsberg, per angular momentum)
and k^shp,l (the shell-polynomial angular-momentum map) -- plus per-element k^diat sigma/pi/
delta frame-scaling constants (Eq 31) and one k^shp element amplitude. The R=6 null sweep
could NOT see overlap-shaped responses (S~ is dead there), so every global is re-probed at
short R where H0 lives. Expected signatures, written before the run:

  k^W_s      moves ss and spz (constant rho vs the level response), never pi
  k^W_p      moves spz/pzpz/pxpx (constant rho), sigma AND pi alike
  k^shp,l    rho GROWS ~linearly with R (the Pi = 1 + k*R polynomial's derivative)
  k^diat,s/p (element slots, L8 quartet candidates: F [2.11, 2.40, 1.07, 0.15]) -- sigma
             slot moves ss/spz/pzpz only; pi slot moves pxpx ONLY; delta slot dead (no d AOs);
             H's L8[1] = 0.0 is consistent with "H cannot pi-bond"
  (already identified from atlas part 1: L1[7] rho is LINEAR in R on both molecules --
   slope 1.28-1.32 on F2, ~0.5 on H2 -- the k^shp ELEMENT amplitude signature.)

F2 at R = 2.0 and 3.8 (two points give the rho slope); H2 at R = 2.5 (sigma-only cross-check).
Level responses for the rho denominators come from the banked atlas/hmetric runs.
Writes data/atlas-globals.json.
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

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
OUT = os.path.join(HERE, "data", "atlas-globals.json")

ELS_F2 = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}
ELS_H2 = {"ss": (0, 1)}
D = 0.02


def probe(sym, R, els, slots):
    out = {}
    for name, key, val in slots:
        params.write_perturbed({key: val + D})
        Fp = fock_recon.fock_ao([(sym, 0, 0, 0), (sym, 0, 0, R / K.BOHR)])["F"]
        params.write_perturbed({key: val - D})
        Fm = fock_recon.fock_ao([(sym, 0, 0, 0), (sym, 0, 0, R / K.BOHR)])["F"]
        out[name] = {el: float((Fp[ix] - Fm[ix]) / (2 * D)) for el, ix in els.items()}
    return out


def main():
    P = params.parse()
    g1, g2 = P["globals"]["g1"], P["globals"]["g2"]
    glb = [(f"G1[{i}]", ("g", 1, i), g1[i]) for i in range(10)] + \
          [(f"G2[{i}]", ("g", 2, i), g2[i]) for i in range(10)]
    l8f = [(f"L8[{i}]F", (9, 8, i), P["element"][9]["l8"][i]) for i in (0, 1, 2, 3, 4, 7)]
    l8h = [(f"L8[{i}]H", (1, 8, i), P["element"][1]["l8"][i]) for i in (0, 1, 2, 3, 4, 7)]
    res = {}
    try:
        for R in (2.0, 3.8):
            res[f"F2@{R}"] = probe("F", R, ELS_F2, glb + l8f)
            print(f"  F2 R={R} done", flush=True)
        res["H2@2.5"] = probe("H", 2.5, ELS_H2, glb + l8h)
        print("  H2 done", flush=True)
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored", flush=True)

    # level denominators from the banked atlas
    atlas = json.load(open(os.path.join(HERE, "data", "offdiag-atlas.json")))["atlas"]
    lev = {}
    for r in atlas["F2"]["recs"]:
        if r["R"] in (2.0, 3.8):
            lev[f"F2@{r['R']}"] = {el: r["resp"]["L2_s"][el] + r["resp"]["L2_p"][el]
                                   for el in ELS_F2}
    for r in atlas["H2"]["recs"]:
        if r["R"] == 2.5:
            lev["H2@2.5"] = {el: r["resp"]["L2_s"][el] for el in ELS_H2}

    print("\nRHO (response/level-response); F2 columns: R=2.0 then R=3.8 per element:")
    print(f"  {'slot':10s}" + "".join(f"{el + '@2':>10s}{el + '@3.8':>10s}" for el in ELS_F2) +
          f"{'H2ss@2.5':>10s}")
    interesting = []
    for name, _k, _v in glb:
        vals = []
        for el in ELS_F2:
            for tag in ("F2@2.0", "F2@3.8"):
                r_ = res[tag][name][el]
                l_ = lev[tag][el]
                vals.append(r_ / l_ if abs(l_) > 1e-8 else float("nan"))
        h_ = res["H2@2.5"][name]["ss"] / lev["H2@2.5"]["ss"]
        vals.append(h_)
        if any(np.isfinite(v) and abs(v) > 0.02 for v in vals):
            interesting.append(name)
            print(f"  {name:10s}" + "".join(f"{v:+10.4f}" for v in vals))
    print("\nL8 quartet (F2, raw slopes; sigma slot -> ss/spz/pzpz, pi slot -> pxpx only):")
    for name, _k, _v in l8f:
        for tag in ("F2@2.0", "F2@3.8"):
            r_ = res[tag][name]
            print(f"  {name:8s} {tag:7s} " + "  ".join(f"{el} {r_[el]:+8.5f}" for el in ELS_F2))
    print("\nL8 quartet (H2 ss):")
    for name, _k, _v in l8h:
        print(f"  {name:8s} ss {res['H2@2.5'][name]['ss']:+8.5f}")
    json.dump({"responses": res, "level_denominators": lev, "movers": interesting},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
