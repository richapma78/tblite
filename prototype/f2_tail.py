"""f2_tail.py -- the F2 long tail (R to 10), the Euler-closure v2 table, and the pzpz object.

The atlas found the structural law: the off-diagonal core is DEGREE-1 HOMOGENEOUS in the
(level, mu) slots -- Euler over (L2, L7) closes H2 ss to +-0.006 Eh and F2 ss/pxpx/spz to a
few percent -- EXCEPT F2's p-sigma element, which keeps a Coulomb-tailed, L4_p-owned object
that no (L2, mu) subtraction touches. This script:

  1. extends the F2 decomposition to R = 7, 8, 10 (gap-gated; stretched RKS F2 is sick at
     some point -- excluded points are listed, not used),
  2. completes the (L2_s, L2_p, L7_s, L7_p) FDs at every stretch point (the atlas covered 4
     of 10), giving the CLOSURE v2 table: object_el(R) = EHT_el - sum slot*dF/dslot,
  3. isolates the pzpz object over the full range and fits tail candidates on R >= 3.8:
        c/R          (bare Coulomb)
        c*erf(aR)/R  (EEQ-style attenuated Coulomb)
        c/(R + b)    (shifted Coulomb / gamma-kernel flavor)
     Pre-declared: a candidate WINS iff its rms <= 2e-3 Eh AND is 2x under the runner-up;
     otherwise the numbers are reported and nothing is claimed.

The object's absolute amplitude carries the known ambiguity of a possible hardcoded
(slot-free) h-bar part, which is S_eff-shaped and so cannot fake the tail. Writes
data/f2-tail.json.
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
import f2_stretch  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
OUT = os.path.join(HERE, "data", "f2-tail.json")

NEW_RS = [7.0, 8.0, 10.0]
ELS = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}
HSLOTS = (("L2_s", 2, 0), ("L2_p", 2, 1), ("L7_s", 7, 0), ("L7_p", 7, 1))
D = 0.02


def hbar_fd(R):
    """(L2/L7 s,p) responses for all elements at R."""
    e = params.parse()["element"][9]
    out = {}
    for name, row, col in HSLOTS:
        v0 = e["shells"][row - 2][col]
        params.write_perturbed({(9, row, col): v0 + D})
        Fp = fock_recon.fock_ao(f2_stretch.f2(R))["F"]
        params.write_perturbed({(9, row, col): v0 - D})
        Fm = fock_recon.fock_ao(f2_stretch.f2(R))["F"]
        out[name] = {el: float((Fp[ix] - Fm[ix]) / (2 * D)) for el, ix in ELS.items()}
    return out


def main():
    e = params.parse()["element"][9]
    slotv = {"L2_s": e["shells"][0][0], "L2_p": e["shells"][0][1],
             "L7_s": e["shells"][5][0], "L7_p": e["shells"][5][1]}
    L5s, L5p = e["shells"][3][0], e["shells"][3][1]
    rows = json.load(open(os.path.join(HERE, "data", "f2-stretch.json")))["rows"]
    atlas = json.load(open(os.path.join(HERE, "data", "offdiag-atlas.json")))
    at_resp = {r["R"]: r["resp"] for r in atlas["atlas"]["F2"]["recs"]}
    tail_rows, resp_all, sick = [], {}, []
    try:
        # ---- 1. new long-R decomposition rows
        for R in NEW_RS:
            params.write_perturbed({})
            r0 = oracle.run(f2_stretch.f2(R))
            if r0["gap_ev"] is None or r0["gap_ev"] < 0.5:
                sick.append({"R": R, "gap_ev": r0["gap_ev"]})
                print(f"  R={R:5.2f} EXCLUDED (gap {r0['gap_ev']})", flush=True)
                continue
            base = fock_recon.fock_ao(f2_stretch.f2(R))
            Fm = {}
            for tag, col, val in (("s+", 0, L5s + 0.05), ("s-", 0, L5s - 0.05),
                                  ("p+", 1, L5p + 0.05), ("p-", 1, L5p - 0.05)):
                params.write_perturbed({(9, 5, col): val})
                Fm[tag] = fock_recon.fock_ao(f2_stretch.f2(R))["F"]
            X = (L5s * (Fm["s+"] - Fm["s-"]) + L5p * (Fm["p+"] - Fm["p-"])) / 0.1
            xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
            A = f2_stretch.acp_matrix([9, 9], xyz)
            EHT = base["F"] - X - A
            row = {"R": R, "gap_ev": r0["gap_ev"],
                   "EHT_el": {el: float(EHT[ix]) for el, ix in ELS.items()}}
            tail_rows.append(row)
            print(f"  R={R:5.2f} gap {r0['gap_ev']:5.2f}  EHT pzpz {row['EHT_el']['pzpz']:+.5f}  "
                  f"ss {row['EHT_el']['ss']:+.5f}", flush=True)
        # ---- 2. hbar FDs everywhere they are missing
        all_rows = rows + tail_rows
        for r in all_rows:
            R = r["R"]
            if R in at_resp:
                resp_all[R] = {n: at_resp[R][n] for n, _row, _col in
                               ((x[0], x[1], x[2]) for x in HSLOTS)}
            else:
                resp_all[R] = hbar_fd(R)
                print(f"  hbar FDs done at R={R:5.2f}", flush=True)
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored", flush=True)

    # ---- closure v2
    print("\nCLOSURE v2: object_el = EHT - Euler(L2_s,L2_p,L7_s,L7_p):")
    print("   R      " + "".join(f"{el:>11s}" for el in ELS) + "    pzpz_obj*R")
    table = []
    for r in sorted(all_rows, key=lambda x: x["R"]):
        R = r["R"]
        rec = {"R": R, "object": {}}
        for el in ELS:
            euler = sum(slotv[n] * resp_all[R][n][el] for n in slotv)
            rec["object"][el] = r["EHT_el"][el] - euler
        table.append(rec)
        print(f"  {R:5.2f} " + "".join(f"{rec['object'][el]:+11.5f}" for el in ELS) +
              f"   {rec['object']['pzpz'] * R:+9.5f}")

    # ---- 3. pzpz tail fits
    pts = [(rec["R"], rec["object"]["pzpz"]) for rec in table if rec["R"] >= 3.8]
    Rv = np.array([p[0] for p in pts])
    Yv = np.array([p[1] for p in pts])

    def fit(fn, grid):
        best = None
        for prm in grid:
            model = fn(Rv, prm)
            c = float(model @ Yv / (model @ model))
            rms = float(np.sqrt(np.mean((c * model - Yv) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, c, prm)
        return best

    cands = {
        "c/R": fit(lambda R, _p: 1 / R, [None]),
        "c*erf(aR)/R": fit(lambda R, a: np.array([math.erf(a * x) for x in R]) / R,
                           np.arange(0.05, 2.01, 0.01)),
        "c/(R+b)": fit(lambda R, b: 1 / (R + b), np.arange(-2.0, 6.01, 0.05)),
    }
    print(f"\npzpz object tail fits (R >= 3.8, {len(pts)} pts):")
    for name, (rms, c, prm) in sorted(cands.items(), key=lambda t: t[1][0]):
        ptxt = "" if prm is None else f"  param {prm:+.2f}"
        print(f"  {name:12s} rms {rms:.2e}  c {c:+.4f}{ptxt}")
    ranked = sorted(cands.items(), key=lambda t: t[1][0])
    win = ranked[0]
    decided = win[1][0] <= 2e-3 and win[1][0] * 2 <= ranked[1][1][0]
    print(f"  -> {'WINNER: ' + win[0] if decided else 'NO WINNER (pre-declared bar not met)'}")

    json.dump({"tail_rows": tail_rows, "sick": sick,
               "closure_v2": table,
               "fits": {n: {"rms": v[0], "c": v[1],
                            "param": (None if v[2] is None else float(v[2]))}
                        for n, v in cands.items()},
               "winner": win[0] if decided else None},
              open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
