"""hmetric.py -- measure the Hamiltonian's OWN metric by level-FD, no basis guessing.

The decoded off-diagonal structure says the level enters as H0_munu ~ K * h-bar * (metric), so
the response of a reconstructed off-diagonal Fock element to the L2 level slot,

    S_eff_munu(R) = -dF_munu/dL2_lbar / K        (K = 2.26, the measured Hueckel prefactor;
                                                  the element's slot shifts BOTH centers of a
                                                  homonuclear pair, exactly as K was measured)

IS the Hamiltonian's effective overlap for that shell pair -- read out of the binary itself.
No scaled-basis scan, no assumed functional form. Then two questions per element type:

  1. SHAPE: is S_eff Gaussian-overlap-like (exponential-ish tail), or does it carry the
     Coulombic tail the F2 p-sigma element shows? Compare against the primary overlap, scaled
     overlaps, and 1/R.
  2. CLOSURE: does EHT_el / S_eff_el go FLAT across the whole stretch (h-bar constant per
     element at fixed shells), including the short range that probe 3 could not explain?

Runs the F2 stretch (L2_s and L2_p FDs) and the H2 stretch (L2_s FD) with the F2/eht_h2
decompositions' banked EHT values read from their JSONs. Writes data/hmetric.json.
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
OUT = os.path.join(HERE, "data", "hmetric.json")
D = 0.02
KHUECKEL = 2.26                       # measured, push 33 (banked h0_offdiagonal_composition)

F2_ELS = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}
H2_ELS = {"ss": (0, 1)}


def diatomic(sym, R_bohr):
    return [(sym, 0, 0, 0), (sym, 0, 0, R_bohr / K.BOHR)]


def level_fd(sym, z, R, slot_col, els):
    """dF_el/dL2_slot by central FD; returns {el: slope}."""
    e = params.parse()["element"][z]
    v0 = e["shells"][0][slot_col]
    params.write_perturbed({(z, 2, slot_col): v0 + D})
    Fp = fock_recon.fock_ao(diatomic(sym, R))["F"]
    params.write_perturbed({(z, 2, slot_col): v0 - D})
    Fm = fock_recon.fock_ao(diatomic(sym, R))["F"]
    return {el: float((Fp[ix] - Fm[ix]) / (2 * D)) for el, ix in els.items()}


def main():
    out = {"K_used": KHUECKEL, "f2": [], "h2": []}
    try:
        f2rows = json.load(open(os.path.join(HERE, "data", "f2-stretch.json")))["rows"]
        for r in f2rows:
            R = r["R"]
            sl_s = level_fd("F", 9, R, 0, F2_ELS)
            sl_p = level_fd("F", 9, R, 1, F2_ELS)
            rec = {"R": R,
                   "Seff_ss": -sl_s["ss"] / KHUECKEL,
                   "Seff_spz_via_s": -sl_s["spz"] / KHUECKEL,
                   "Seff_spz_via_p": -sl_p["spz"] / KHUECKEL,
                   "Seff_pzpz": -sl_p["pzpz"] / KHUECKEL,
                   "Seff_pxpx": -sl_p["pxpx"] / KHUECKEL,
                   "cross_s_on_pzpz": -sl_s["pzpz"] / KHUECKEL,
                   "EHT": r["EHT_el"], "S_primary": {"ss": r["S_ss"], "spz": r["S_spz"],
                                                     "pzpz": r["S_pzpz"], "pxpx": r["S_pxpx"]}}
            out["f2"].append(rec)
            print(f"  F2 R={R:4.2f}  Seff ss {rec['Seff_ss']:+.5f}  "
                  f"spz {rec['Seff_spz_via_s']:+.5f}/{rec['Seff_spz_via_p']:+.5f}  "
                  f"pzpz {rec['Seff_pzpz']:+.5f}  pxpx {rec['Seff_pxpx']:+.5f}", flush=True)
        h2rows = json.load(open(os.path.join(HERE, "data", "h2-eht-probes.json")))["rows"]
        for r in h2rows:
            R = r["R"]
            sl = level_fd("H", 1, R, 0, H2_ELS)
            rec = {"R": R, "Seff_ss": -sl["ss"] / KHUECKEL, "EHT_ss": r["EHT12"],
                   "S_primary": r["S12"]}
            out["h2"].append(rec)
            print(f"  H2 R={R:4.2f}  Seff {rec['Seff_ss']:+.5f}  S {r['S12']:+.5f}  "
                  f"EHT/Seff {r['EHT12'] / rec['Seff_ss'] if rec['Seff_ss'] else float('nan'):+.4f}",
                  flush=True)
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored", flush=True)

    # ---------------------------------------------------------------- closure + shape report
    print("\nCLOSURE: EHT_el / Seff_el per R (flat = the metric explains the element):")
    hdr = "   R    " + "".join(f"{el:>12s}" for el in ("ss", "spz", "pzpz", "pxpx"))
    print("  F2:" + hdr)
    for rec in out["f2"]:
        vals = []
        for el, seff_key in (("ss", "Seff_ss"), ("spz", "Seff_spz_via_p"),
                             ("pzpz", "Seff_pzpz"), ("pxpx", "Seff_pxpx")):
            s = rec[seff_key]
            vals.append(rec["EHT"][el] / s if abs(s) > 1e-7 else float("nan"))
        print(f"  {rec['R']:5.2f} " + "".join(f"{v:+12.4f}" for v in vals))
    print("  H2:")
    for rec in out["h2"]:
        s = rec["Seff_ss"]
        print(f"  {rec['R']:5.2f}  EHT/Seff {rec['EHT_ss'] / s if abs(s) > 1e-7 else float('nan'):+10.4f}"
              f"   Seff/S_primary {s / rec['S_primary']:+8.4f}")
    print("\nSHAPE (F2 pzpz): Seff vs 1/R vs primary S:")
    for a, b in zip(out["f2"], out["f2"][1:]):
        if abs(a["Seff_pzpz"]) > 1e-7 and abs(b["Seff_pzpz"]) > 1e-7:
            print(f"  {a['R']:4.2f}->{b['R']:4.2f}: Seff ratio {a['Seff_pzpz'] / b['Seff_pzpz']:6.3f}"
                  f"   1/R ratio {b['R'] / a['R']:6.3f}"
                  f"   S ratio {a['S_primary']['pzpz'] / b['S_primary']['pzpz'] if b['S_primary']['pzpz'] else float('nan'):6.3f}")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
