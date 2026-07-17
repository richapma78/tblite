"""offdiag_atlas.py -- the slot-response atlas for off-diagonal Fock elements.

The original perturbation map cracked the parameter file by asking which printed observable
each slot moves. This extends it to the off-diagonal Fock: for each h-bar-candidate slot,
the response dF_el/dslot at several R, SHAPE-NORMALIZED by the same element's level response

    rho_slot(R) = (dF_el/dslot) / (dF_el/dL2bar)          (K and the metric cancel)

Selection rule (pre-declared): a slot with |rho| variation < 25% across the R set and a
response above noise is an h-bar MEMBER with relative weight rho; a slot whose rho GROWS with
R belongs to the second (L4/L1[7]) channel. L4 and L1[7] are tabulated for the record but are
channel knobs by construction (hmetric: the level channel rides the PRIMARY overlap, which no
L4/L1[7] touches -- dF11/dL4 = dF11/dL1[7] = 0 on H2).

Then the EULER SUM over the h-bar family isolates the channel with no model and no fit:

    channel_el(R) = EHT_el(R) - sum_{hbar slots} slot * dF_el/dslot

exact up to a possible HARDCODED h-bar part (S_eff-shaped, constant amplitude). An S_eff-shaped
ambiguity cannot fake a Coulomb tail, so the channel's TAIL FORM is unambiguous; its absolute
amplitude carries the ambiguity and is labeled so.

EHT values come from the banked runs (data/f2-stretch.json, data/h2-eht-probes.json); the
atlas R values are chosen inside those grids. Writes data/offdiag-atlas.json.
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
OUT = os.path.join(HERE, "data", "offdiag-atlas.json")

CASES = {
    "F2": {"sym": "F", "z": 9, "rs": [2.0, 2.668, 3.8, 5.2],
           "els": {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)},
           "banked": "f2-stretch.json", "shells": ("s", "p")},
    "H2": {"sym": "H", "z": 1, "rs": [1.4, 2.5, 4.0, 6.0],
           "els": {"ss": (0, 1)}, "banked": "h2-eht-probes.json", "shells": ("s",)},
}


def slot_list(z, shells):
    e = params.parse()["element"][z]
    out = []
    for li, (row, name) in enumerate(((2, "L2"), (3, "L3"), (7, "L7"))):
        for ci, sh in enumerate(shells):
            out.append((f"{name}_{sh}", (z, row, ci), e["shells"][row - 2][ci], 0.02))
    out.append(("L1[8]k1cn", (z, 1, 8), e["l1"][8], 0.02))
    out.append(("L1[7]hbas", (z, 1, 7), e["l1"][7], 0.01))
    for ci, sh in enumerate(shells):
        out.append((f"L4_{sh}", (z, 4, ci), e["shells"][2][ci], 0.05))
    return out


def diatomic(sym, R):
    return [(sym, 0, 0, 0), (sym, 0, 0, R / K.BOHR)]


def main():
    atlas = {}
    try:
        for case, cfg in CASES.items():
            slots = slot_list(cfg["z"], cfg["shells"])
            recs = []
            for R in cfg["rs"]:
                rec = {"R": R, "resp": {}}
                for name, key, val, d in slots:
                    params.write_perturbed({key: val + d})
                    Fp = fock_recon.fock_ao(diatomic(cfg["sym"], R))["F"]
                    params.write_perturbed({key: val - d})
                    Fm = fock_recon.fock_ao(diatomic(cfg["sym"], R))["F"]
                    rec["resp"][name] = {el: float((Fp[ix] - Fm[ix]) / (2 * d))
                                         for el, ix in cfg["els"].items()}
                recs.append(rec)
                lv = {el: rec["resp"].get("L2_s", {}).get(el, 0) +
                      rec["resp"].get("L2_p", {}).get(el, 0) for el in cfg["els"]}
                print(f"  [{case}] R={R:4.2f} done (level resp ss {lv['ss']:+.4f})", flush=True)
            atlas[case] = {"slots": [(n, k, v, d) for n, k, v, d in slots], "recs": recs}
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored", flush=True)

    # ------------------------------------------------------------------- rho tables
    print("\nRHO TABLES: rho = resp / (resp_L2_s + resp_L2_p) per element; constant = hbar "
          "member, growing = channel:")
    hbar_family = {}
    for case, cfg in CASES.items():
        recs = atlas[case]["recs"]
        for el in cfg["els"]:
            print(f"  [{case}:{el}]")
            lead = [sum(r["resp"][f"L2_{sh}"][el] for sh in cfg["shells"]) for r in recs]
            fam = []
            for name, _k, val, _d in atlas[case]["slots"]:
                rhos = [r["resp"][name][el] / l if abs(l) > 1e-8 else float("nan")
                        for r, l in zip(recs, lead)]
                a = [abs(x) for x in rhos if np.isfinite(x)]
                var = (max(a) - min(a)) / max(a) if a and max(a) > 1e-4 else None
                is_ch = name.startswith(("L4", "L1[7]"))
                tag = ("channel-knob" if is_ch else
                       "hbar" if var is not None and var < 0.25 else
                       "grows/odd" if var is not None else "quiet")
                if tag == "hbar":
                    fam.append(name)
                print(f"    {name:10s} rho " + " ".join(f"{x:+8.4f}" for x in rhos) +
                      f"   [{tag}]")
            hbar_family[(case, el)] = fam

    # ------------------------------------------------------------ channel extraction
    print("\nCHANNEL: EHT - Euler(hbar family incl. levels) per element "
          "(+ unknown S_eff-shaped constant):")
    banked = {c: json.load(open(os.path.join(HERE, "data", CASES[c]["banked"])))["rows"]
              for c in CASES}
    chan_out = {}
    for case, cfg in CASES.items():
        recs = atlas[case]["recs"]
        slots = {n: (v) for n, _k, v, _d in atlas[case]["slots"]}
        for el in cfg["els"]:
            fam = [f"L2_{sh}" for sh in cfg["shells"]] + hbar_family[(case, el)]
            fam = list(dict.fromkeys(fam))
            rowsB = banked[case]
            line = []
            for rec in recs:
                bank = min(rowsB, key=lambda r: abs(r["R"] - rec["R"]))
                eht = bank["EHT_el"][el] if case == "F2" else bank["EHT12"]
                euler = sum(slots[n] * rec["resp"][n][el] for n in fam)
                ch = eht - euler
                line.append({"R": rec["R"], "EHT": eht, "euler_hbar": euler, "channel": ch})
            chan_out[f"{case}:{el}"] = {"family": fam, "points": line}
            print(f"  [{case}:{el}] family {fam}")
            for p in line:
                print(f"    R={p['R']:4.2f}  EHT {p['EHT']:+.5f}  hbar-part {p['euler_hbar']:+.5f}"
                      f"  channel {p['channel']:+.5f}  channel*R {p['channel'] * p['R']:+.5f}")
    json.dump({"atlas": atlas, "channel": chan_out}, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
