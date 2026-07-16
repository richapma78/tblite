"""chart_onsite.py -- the dense charting campaign for the HARDCODED onsite functions.

Three unknown functions block every molecular gate (proven hardcoded in the binary):
  s(...)  -- the scale on the harmonic-kernel onsite 2nd order
  E3 rule -- the cubic term's dependence on the Hubbards (Gamma element rule)
  E4(...) -- the quartic-ish onsite remainder

The instrument (SET, don't nudge): for each charge state of an element, SET the L6 Hubbards to
chosen values and read the printed ES2+3.
  - equal-U scans at several u separate E2 (known shape, ~q_A^2*u at equal U... exactly
    1/2*q_A^2*s*u) from E3 (u-dependence reveals whether Gamma scales with U) and E4.
  - asymmetry scans (U_s grid, U_p fixed) give the harmonic-shape coefficient K = s * f(q)
    per charge state -> s as a function of the charges.
Shell charges per state come from the printed populations + the fractional references
(constants.REFOCC), aufbau-filled ions verified by the FD technique earlier.

Writes data/onsite-charts.json (raw scans + derived quantities). Analysis prints inline.
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

CASES = {                       # element -> {chg: uhf}
    "o": {"z": 8, "states": {-2: 0, -1: 1, 0: 2, 1: 3, 2: 2}},
    "c": {"z": 6, "states": {-2: 0, -1: 1, 0: 2, 1: 1, 2: 0}},
}
EQUAL_US = (0.5, 1.0, 2.0)
ASYM_US = (0.5, 1.0, 1.5, 2.0, 3.0)


def _lines():
    return open(PRISTINE).read().splitlines(keepends=True)


def _l6_line(lines, z):
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == [str(z)])
    return idx[idx.index(hdr) + 6]


def set_l6(lines, z, us, up):
    li = _l6_line(lines, z)
    t = lines[li].split()
    t[0], t[1] = f"{us:.10f}", f"{up:.10f}"
    out = list(lines)
    out[li] = "      " + "      ".join(t) + "\n"
    open(HOME, "w").writelines(out)


def run(sym, chg, uhf):
    r = oracle.run([(sym, 0.0, 0.0, 0.0)], charge=chg, uhf=uhf)
    pops = r["pops"][0]
    return r["terms"]["ES2+3"], (pops["s"], pops["p"])


def main():
    lines = _lines()
    data = {}
    try:
        for sym, spec in CASES.items():
            z = spec["z"]
            ref = K.REFOCC[z]
            for chg, uhf in spec["states"].items():
                key = f"{sym}{chg:+d}"
                rec = {"z": z, "chg": chg, "equal": {}, "asym": {}}
                shutil.copyfile(PRISTINE, HOME)
                es0, (ps, pp) = run(sym, chg, uhf)
                qs, qp = ref[0] - ps, ref[1] - pp
                rec["q"] = [qs, qp]
                rec["pristine_es23"] = es0
                for u in EQUAL_US:
                    set_l6(lines, z, u, u)
                    e, _ = run(sym, chg, uhf)
                    rec["equal"][str(u)] = e
                for us in ASYM_US:
                    set_l6(lines, z, us, 1.0)
                    e, _ = run(sym, chg, uhf)
                    rec["asym"][str(us)] = e
                data[key] = rec
                print(f"  {key}: q=({qs:+.4f},{qp:+.4f})  equal-U " +
                      " ".join(f"{u}:{rec['equal'][str(u)]:+.6f}" for u in EQUAL_US), flush=True)
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("restored")
    json.dump(data, open(os.path.join(HERE, "data", "onsite-charts.json"), "w"), indent=1)

    # ------------------------------------------------------------------ inline analysis
    print("\nANALYSIS")
    print("1) E3 u-dependence (odd part of equal-U ion pairs, per u):")
    for sym in CASES:
        for u in EQUAL_US:
            su = str(u)
            for a in (1, 2):
                kp, km = f"{sym}{a:+d}", f"{sym}{-a:+d}"
                if kp in data and km in data:
                    odd = 0.5 * (data[kp]["equal"][su] - data[km]["equal"][su])
                    print(f"   {sym} |chg|={a} u={u}: odd = {odd:+.6f}  odd/q^3 = "
                          f"{odd / a**3:+.6f}")
    print("2) asymmetry K (fit E = C + K*(us-up)^2/(us+up), up=1) -> per state:")
    for key, rec in data.items():
        xs = [(us, rec["asym"][str(us)]) for us in ASYM_US]
        # linear fit E vs t=(us-1)^2/(us+1)
        ts = np.array([(us - 1.0) ** 2 / (us + 1.0) for us, _ in xs])
        es = np.array([e for _, e in xs])
        A = np.vstack([ts, np.ones_like(ts)]).T
        (Kfit, Cfit), res, *_ = np.linalg.lstsq(A, es, rcond=None)
        resid = float(np.abs(A @ np.array([Kfit, Cfit]) - es).max())
        qs, qp = rec["q"]
        # with q_p relative... harmonic-shape coefficient for general (qs, qp):
        # E2 = s*[1/2 qs^2 Us + 1/2 qp^2 Up + qs qp gam_harm]; against t-shape only the
        # (qs qp) cross term carries the (dU)^2/(sumU) form when Us,Up vary... report raw K.
        print(f"   {key}: K = {Kfit:+.6f}  C = {Cfit:+.6f}  shape-resid {resid:.1e}  "
              f"q=({qs:+.3f},{qp:+.3f})  K/(-qs*qp) = "
              f"{Kfit / (-qs * qp) if abs(qs * qp) > 1e-9 else float('nan'):+.4f}")


if __name__ == "__main__":
    main()
