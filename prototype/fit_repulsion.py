"""fit_repulsion.py -- identify the repulsion parameters (SI Sec. 1.6) in gxtb_parameters.

Model (SI Eqs. 52-56), homonuclear neutral pair (q = 0 -> Zeff = Zeff0; alpha_AB = alpha_A/2):

  E(R) = Zeff0^2 * (1 + p1/R + p2/R^2 + p3/R^3 + p4/R^4) * exp(-aeff * (R (+/-) R0)^1.5)
  aeff = (alpha0/2) * (1 + kCN*sqrt(CN))-ish -- CN ~ 0 in the tail, so tail fits (alpha0/2, R0).

Hypotheses under test (from the cross-element column trends):
  Zeff0 = line1 col3 (grows with Z: H 1.17 ... I 7.8); alpha0, R0 among cols 1/6;
  p1..p4 among the 20 global floats. Strategy: fix Zeff0 = col3, fit the H2 tail (R >= 3,
  CN ~ 0) for (p1, aeff, R0) under both offset signs, then compare fitted values against the
  file's candidate columns and the globals -- the find-the-float trick, again.
"""
import itertools
import math
import sys

import numpy as np

# oracle H2 nuclear-repulsion curve (R_bohr -> Eh), measured 2026-07-16
CURVE = {1.0: 0.27544985, 1.2: 0.17143088, 1.4: 0.10978789, 1.6: 0.06983107,
         2.0: 0.02661631, 2.5: 0.00702091, 3.0: 0.00161783, 4.0: 0.00005958,
         5.0: 0.00000142, 6.0: 0.00000002}
ZEFF_H = 1.1699036268          # line1 col3 hypothesis
H_LINE1 = [1.1527189408, 0.2105955739, 1.1699036268, -0.0067641524, -0.2141751345,
           1.0546334924, 0.8141717488, 0.0223676961, 0.7749747393, -0.1978436058]
GLOBALS = [0.7188463521, 0.8284173386, 1.19953196, 0.2599875891, 0.0788775224,
           1.7995847408, -0.5959929766, 0.2140456651, -0.3191336623, 1.2154627292,
           -1.8217544049, 0.3300126723, -0.0276743359, 0.1057773939, 0.1682969331,
           0.4667106684, 0.2287119112, 0.0522538831, 0.2347047181, 0.304294728]


def model(R, p1, aeff, r0, sign, p2=0.0, p3=0.0, p4=0.0):
    pen = 1.0 + p1 / R + p2 / R**2 + p3 / R**3 + p4 / R**4
    return ZEFF_H**2 * pen * math.exp(-aeff * (R + sign * r0) ** 1.5)


def fit_tail(sign):
    """Grid + refine on (p1, aeff, r0) over the tail R >= 2.5 (log residuals)."""
    Rs = [2.5, 3.0, 4.0, 5.0, 6.0]
    ys = [math.log(CURVE[r]) for r in Rs]

    def resid(p1, aeff, r0):
        try:
            return sum((math.log(model(r, p1, aeff, r0, sign)) - y) ** 2
                       for r, y in zip(Rs, ys))
        except ValueError:
            return 1e9

    best = None
    for p1 in np.arange(-2.0, 2.01, 0.25):
        for aeff in np.arange(0.2, 1.51, 0.05):
            for r0 in np.arange(0.0, 2.01, 0.1):
                r = resid(p1, aeff, r0)
                if best is None or r < best[0]:
                    best = (r, p1, aeff, r0)
    # local refine
    _, p1, aeff, r0 = best
    step = 0.05
    for _ in range(60):
        improved = False
        for dp, da, dr in itertools.product((-step, 0, step), repeat=3):
            r = resid(p1 + dp, aeff + da, r0 + dr)
            if r < best[0]:
                best = (r, p1 + dp, aeff + da, r0 + dr)
                improved = True
        _, p1, aeff, r0 = best
        if not improved:
            step /= 2
            if step < 1e-7:
                break
    return best


def main():
    for sign, name in ((+1, "exp(-a(R+R0)^1.5)"), (-1, "exp(-a(R-R0)^1.5)")):
        res, p1, aeff, r0 = fit_tail(sign)
        print(f"  {name}:  res {res:.2e}   p1 {p1:+.4f}  aeff {aeff:.4f} "
              f"(alpha0 = {2*aeff:.4f})  R0 {r0:.4f}")
        # nearest candidates in the file
        cand_a = min(H_LINE1, key=lambda v: abs(v - 2 * aeff))
        cand_r = min(H_LINE1, key=lambda v: abs(v - r0))
        cand_p = min(GLOBALS, key=lambda v: abs(v - p1))
        print(f"      nearest file floats: alpha0~{cand_a}  R0~{cand_r}  p1(global)~{cand_p}")


if __name__ == "__main__":
    sys.exit(main())
