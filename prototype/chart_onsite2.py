"""chart_onsite2.py -- campaign v2: per-point POPULATIONS, s-rule subtracted, E3 shape isolated.

v1's lesson: ion scans shuffle electrons between shells as the Hubbards change, so any analysis
assuming fixed shell charges is contaminated (Ge's wild s, F's k3 anomaly). Fix: read the printed
populations AT EVERY SCAN POINT -- q_l(point) is then exact -- and subtract the now-KNOWN
second-order (s-rule x harmonic kernel) pointwise. The residual is pure third-plus-higher order:

    resid(U_s, U_p; q_l) = ES23 - s*[1/2 sum q_l q_l' gamma_harm(U)]

charted over the asymmetry grid for each ion state. Writes data/onsite-charts2.json;
analysis prints the residual structure per state.
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

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")

S_RULE = {"period2": (0.08247, 0.0920), "period3": (0.04193, 0.1378)}
NVAL = {6: 4, 7: 5, 8: 6, 9: 7, 14: 4, 15: 5, 16: 6, 17: 7}
PERIOD = {6: 2, 7: 2, 8: 2, 9: 2, 14: 3, 15: 3, 16: 3, 17: 3}

CASES = {
    "o": (8, {-1: 1, 1: 3, 2: 2}),
    "c": (6, {-1: 1, 1: 1}),
    "f": (9, {-1: 0, 1: 2}),
    "s": (16, {-1: 1, 1: 3}),
}
GRID = [(0.5, 0.5), (1.0, 1.0), (2.0, 2.0), (0.5, 1.0), (1.5, 1.0), (2.0, 1.0), (3.0, 1.0),
        (1.0, 0.5), (1.0, 2.0)]


def s_of(z):
    a, b = S_RULE[f"period{PERIOD[z]}"]
    return a * NVAL[z] + b


def _lines():
    return open(PRISTINE).read().splitlines(keepends=True)


def set_l6(lines, z, us, up):
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == [str(z)])
    li = idx[idx.index(hdr) + 6]
    t = lines[li].split()
    t[0], t[1] = f"{us:.10f}", f"{up:.10f}"
    out = list(lines)
    out[li] = "      " + "      ".join(t) + "\n"
    open(HOME, "w").writelines(out)


def main():
    lines = _lines()
    data = {}
    try:
        for sym, (z, states) in CASES.items():
            ref = K.REFOCC[z]
            s = s_of(z)
            for chg, uhf in states.items():
                key = f"{sym}{chg:+d}"
                pts = []
                for us, up in GRID:
                    set_l6(lines, z, us, up)
                    r = oracle.run([(sym, 0.0, 0.0, 0.0)], charge=chg, uhf=uhf)
                    p = r["pops"][0]
                    qs, qp = ref[0] - p["s"], ref[1] - p["p"]
                    qd = ref.get(2, 0.0) - p["d"]
                    es = r["terms"]["ES2+3"]
                    q = [qs, qp] + ([qd] if 2 in ref else [])
                    U = [us, up] + ([us] if 2 in ref else [])
                    E2 = s * sum(0.5 * q[a] * q[b] * (2 * U[a] * U[b] / (U[a] + U[b]))
                                 for a in range(len(q)) for b in range(len(q)))
                    pts.append({"us": us, "up": up, "q": q, "es23": es,
                                "resid": es - E2})
                data[key] = {"z": z, "chg": chg, "s_used": s, "points": pts}
                qa = sum(pts[1]["q"])
                print(f"  {key}: q_A ~ {qa:+.3f}  resid at (1,1)/(2,2)/(3,1): " +
                      "  ".join(f"{p['resid']:+.6f}" for p in pts if
                                (p['us'], p['up']) in ((1.0, 1.0), (2.0, 2.0), (3.0, 1.0))),
                      flush=True)
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("restored")
    json.dump(data, open(os.path.join(HERE, "data", "onsite-charts2.json"), "w"), indent=1)

    print("\nresidual structure per state (resid = ES23 - s*E2_harm, pointwise q):")
    for key, rec in data.items():
        qa3 = sum(rec["points"][1]["q"]) ** 3
        print(f"  {key}:")
        for p in rec["points"]:
            us, up = p["us"], p["up"]
            uterm = 0.25 * (us + up) ** 2          # the offsite Eq.132 (U+U')^2/4 candidate
            print(f"    ({us:3.1f},{up:3.1f}) resid {p['resid']:+9.6f}   "
                  f"resid/q_A^3 {p['resid'] / qa3 if abs(qa3) > 1e-9 else float('nan'):+9.6f}   "
                  f"/(q_A^3*(U+U')^2/4) "
                  f"{p['resid'] / (qa3 * uterm) if abs(qa3 * uterm) > 1e-12 else float('nan'):+9.6f}")


if __name__ == "__main__":
    main()
