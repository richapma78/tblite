"""d4_probe.py -- stock-dftd4 vs the binary's printed revD4 dispersion (gap #6).

Runs under ~/dftd4env/bin/python3 (dftd4-python 4.2.0; NOT installed in gpudft).

DERIVED (2026-07-18, falsifiable): the param-file dispersion slots are
  a1 = g1[9] = 1.2154627292   (BJ damping scale; a2 is ELIMINATED in revD4 -- SI Eq 170,
                               R0 = sqrt(C8/C6), so a1 alone carries the radius)
  s8 = g2[9] = 0.304294728    (s6 = 1 presumed, s9 ~ irrelevant here: ATM ~0.003 mEh)
Evidence: stock d4 under this assignment lands within 0.30 mEh of the printed value on all
six references; the REVERSED assignment misses by -100..-2540 mEh.

Stock does NOT close revD4 (d4 0.15-0.30 mEh, d4s 0.10-0.51, mixed signs): the remainder
is revD4's own charge model -- the sigmoidal zeta (SI Eq 166-169) driven by the CONVERGED
MULLIKEN charges (not EEQ; HF Mulliken +-0.412 vs EEQ +-0.117 is a large zeta shift).
Declared wire bar 0.01 mEh NOT met -> the ledger keeps dispersion MISSING. The close:
implement Eq 166-169 over the stock CP machinery; get zeta point-by-point by INJECTING
synthetic charges into the binary's D4 call under gdb (see derived-constants.json,
"dispersion_revd4_slots_DERIVED").
"""
import math

import numpy as np
from dftd4.interface import DampingParam, DispersionModel

BOHR = 1.8897261254578281
ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
W = [[0.0, 0.0, 0.0],
     [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
     [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
a4 = 1.087 * BOHR / math.sqrt(3)
r_nh = 1.012 * BOHR
st, ct = 0.9262, -0.3770
# printed dispersion, read from the same runs the ledger uses (2026-07-18 audit)
SYSTEMS = {
    "h2": ([1, 1], [[0, 0, 0], [0, 0, 1.4]], -0.000676),
    "f2": ([9, 9], [[0, 0, 0], [0, 0, 2.668]], -0.000882),
    "hf": ([1, 9], [[0, 0, 0], [0, 0, 1.733]], -0.000239),
    "h2o": ([8, 1, 1], W, -0.000776),
    "ch4": ([6, 1, 1, 1, 1],
            [[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4], [-a4, a4, -a4], [-a4, -a4, a4]],
            -0.003457),
    "nh3": ([7, 1, 1, 1],
            [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                            r_nh * st * math.sin(2 * math.pi * k / 3),
                            r_nh * ct] for k in range(3)], -0.002118),
}
A1, S8 = 1.2154627292, 0.304294728


def run(model="d4", s9=1.0, verbose=True):
    worst = 0.0
    for name, (zs, xyz, printed) in SYSTEMS.items():
        m = DispersionModel(np.array(zs), np.array(xyz, float), model=model)
        e = m.get_dispersion(
            DampingParam(s6=1.0, s8=S8, a1=A1, a2=0.0, s9=s9), grad=False)["energy"]
        d = (e - printed) * 1000
        worst = max(worst, abs(d))
        if verbose:
            print(f"  {name:4s} ours {e:+.6f}  printed {printed:+.6f}  d {d:+8.3f} mEh")
    if verbose:
        print(f"  worst |d| = {worst:.3f} mEh  (stock-{model} distance to revD4; "
              f"wire bar 0.01 -- NOT met, by design: the zeta model differs)")
    return worst


if __name__ == "__main__":
    run()
