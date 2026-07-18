"""mfx_gamma_gate.py -- regression gate: mfx.gamma_matrix must reproduce the binary's own
AO exchange gamma (4th arg of setgab_lrao_) bit-for-bit across an H2 bond sweep.

The reference p4(R) values were gdb-extracted from /opt/gxtb-v1/gxtb (setgab_lrao_ 0x539ee0,
save $rcx, finish, read the matrix). This gate needs NO binary -- it pins the generation so a
future edit to the kernel or constants that breaks the match fails loudly. Onsite is R-free
(a pure atomic constant); offsite rises-then-falls (range-separated erf numerator).
"""
import numpy as np

import mfx

# (R_bohr, onsite p4[0,0], offsite p4[0,1]) extracted from the binary
REF = [
    (1.0, 0.3601297506, 0.0396855562),
    (1.4, 0.3601297506, 0.0445643577),
    (1.8, 0.3601297506, 0.0481466063),
    (2.4, 0.3601297506, 0.0514631945),
    (3.0, 0.3601297506, 0.0527750600),
    (4.0, 0.3601297506, 0.0517689203),
    (6.0, 0.3601297506, 0.0438507394),
]
# H2 AO meta: two s shells (atom, Z, l, _), one per H
META = [(0, 1, 0, None), (1, 1, 0, None)]


def run():
    worst = 0.0
    for R, on_ref, off_ref in REF:
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        gam = mfx.gamma_matrix([1, 1], xyz, META)
        don, doff = abs(gam[0, 0] - on_ref), abs(gam[0, 1] - off_ref)
        worst = max(worst, don, doff)
        print(f"  R={R:<4} onsite {gam[0,0]:.10f} (d {gam[0,0]-on_ref:+.1e})  "
              f"offsite {gam[0,1]:.10f} (d {gam[0,1]-off_ref:+.1e})")
    ok = worst < 1e-9
    print(f"\n  worst |generated - binary p4| = {worst:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() else 1)
