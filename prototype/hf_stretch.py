"""hf_stretch.py -- the HF stretch: the heteronuclear collection instrument.

HF (nsao 5, fully invertible; polar, q ~ +-0.3) is where three things live that the
homonuclear program could not reach:

  1. the heteronuclear AVERAGING WEIGHTS: Eq. 64 says (H_lA + H_lB)/2 -- arithmetic, equal
     weights across two DIFFERENT elements. The L2 FDs on (Hs,Fs) and (Hs,Fpz) measure the
     weights directly (on F2's spz they summed to 1.0000 with a 0.48/0.52 split).
  2. the POLAR channels: at q != 0 the Ham basis's k-tilde-0 (q) and k-tilde-3 (q*CN)
     adaptations wake up, the f(q) factors on mu turn on, and the ES2 off-diagonal
     (Mulliken 1/2*S*(gamma*q sums)) is nonzero for the first time.
  3. the RESIDUAL AS MEASUREMENT: EHT - (level part) - (mu part) - (ACP/X, subtracted as
     always) = the polar-ES off-diagonal + any q-driven object -- the direct measurement of
     the half-decoded offsite gamma kernel on the Fock.

This script COLLECTS (oracle FDs; the analysis is a separate offline script):
  per R: base reconstruction; X-Euler slots L5(H,s)/L5(F,s)/L5(F,p); level slots
  L2(H,s)/L2(F,s)/L2(F,p); mu slots L7(H,s)/L7(F,s)/L7(F,p); Hubbard slots L6(F,s)/L6(F,p)
  (the ES2 knobs -- gamma is U-dependent); ACP analytic. Elements tracked: (Hs,Fs) = (0,1),
  (Hs,Fpz) = (0,4), plus the onsite-F (Fs,Fpz) = (1,4) null and the pi zeros as checks.
Writes data/hf-stretch.json. Gap-gated; the EEQ charges per point are recorded (the q(R)
curve matters for every polar channel).
"""
import json
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

PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
OUT = os.path.join(HERE, "data", "hf-stretch.json")
RS = [1.4, 1.6, 1.733, 2.0, 2.3, 2.7, 3.2, 3.8, 4.5]
ELS = {"HsFs": (0, 1), "HsFpz": (0, 4), "onsite_FsFpz": (1, 4)}
D = 0.02
FD_SLOTS = [("L5_Hs", (1, 5, 0), 0.05), ("L5_Fs", (9, 5, 0), 0.05), ("L5_Fp", (9, 5, 1), 0.05),
            ("L2_Hs", (1, 2, 0), D), ("L2_Fs", (9, 2, 0), D), ("L2_Fp", (9, 2, 1), D),
            ("L7_Hs", (1, 7, 0), D), ("L7_Fs", (9, 7, 0), D), ("L7_Fp", (9, 7, 1), D),
            ("L6_Fs", (9, 6, 0), 0.05), ("L6_Fp", (9, 6, 1), 0.05)]


def hf(R):
    return [("H", 0, 0, 0), ("F", 0, 0, R / K.BOHR)]


def main():
    P = params.parse()
    vals = {}
    for name, (z, row, col), d in FD_SLOTS:
        e = P["element"][z]
        vals[name] = (e["l1"][col] if row == 1 else
                      e["shells"][row - 2][col] if 2 <= row <= 7 else
                      e["l8"][col] if row == 8 else e["l9"][col])
    L5v = {"L5_Hs": vals["L5_Hs"], "L5_Fs": vals["L5_Fs"], "L5_Fp": vals["L5_Fp"]}
    rows, sick = [], []
    try:
        for R in RS:
            params.write_perturbed({})
            r0 = oracle.run(hf(R))
            if r0["gap_ev"] is None or r0["gap_ev"] < 0.5:
                sick.append({"R": R, "gap_ev": r0["gap_ev"]})
                print(f"  R={R:5.2f} EXCLUDED (gap {r0['gap_ev']})", flush=True)
                continue
            base = fock_recon.fock_ao(hf(R))
            resp = {}
            for name, key, d in FD_SLOTS:
                params.write_perturbed({key: vals[name] + d})
                Fp = fock_recon.fock_ao(hf(R))["F"]
                params.write_perturbed({key: vals[name] - d})
                Fm = fock_recon.fock_ao(hf(R))["F"]
                resp[name] = {el: float((Fp[ix] - Fm[ix]) / (2 * d))
                              for el, ix in ELS.items()}
            X = {el: sum(L5v[s] * resp[s][el] for s in L5v) for el in ELS}
            xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
            A = f2_stretch.acp_matrix([1, 9], xyz)
            row = {"R": R, "gap_ev": r0["gap_ev"],
                   "q_H": r0["eeq"][0]["q"], "q_F": r0["eeq"][1]["q"],
                   "cn_H": r0["eeq"][0]["cn"],
                   "S": {el: float(base["S"][ix]) for el, ix in ELS.items()},
                   "F": {el: float(base["F"][ix]) for el, ix in ELS.items()},
                   "A": {el: float(A[ix]) for el, ix in ELS.items()},
                   "X": X,
                   "EHT": {el: float(base["F"][ix]) - X[el] - float(A[ix])
                           for el, ix in ELS.items()},
                   "resp": resp,
                   "pi_leak": float(max(abs(base["F"][0, 2]), abs(base["F"][0, 3])))}
            rows.append(row)
            print(f"  R={R:5.2f} gap {r0['gap_ev']:5.2f}  q_H {row['q_H']:+.3f}  "
                  f"EHT HsFs {row['EHT']['HsFs']:+.5f}  HsFpz {row['EHT']['HsFpz']:+.5f}  "
                  f"pi-leak {row['pi_leak']:.1e}", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored", flush=True)
    json.dump({"rows": rows, "sick": sick}, open(OUT, "w"), indent=1)
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
