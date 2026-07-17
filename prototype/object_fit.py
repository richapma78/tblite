"""object_fit.py -- the pzpz object closed form (sign-honest fits, reproducible).

RESULT (banked): obj(R) = c*erf(a*R)/R - A*S~(0.635) with c = 0.42, a = 0.23 (= G2[8] =
0.2347, FD-shape corroborated at R=2.67/6), A = 0.048; rms 4.2e-4 over the positive branch
(R = 2..6) and 4.9e-4 on |obj| over all 13 points -- both under the pre-declared 2e-3 bar.
The R>=7 sign flip is an ATTRIBUTION artifact: the crossing-capable fit chose no crossing.
c ~ 3*|drho_p(F)| = 0.432 (3%) is FLAGGED; it survives its first falsification (REFOCC(H) =
1.0 exactly -> drho(H) = 0 -> the H2-absence of the object is predicted). The 6->7 jump (+0.064 -> -0.057 with smooth magnitude)
is not a smooth crossing; the R>=7 SIGN is treated as unreliable (near-degenerate sigma pair;
the R=6 pairing check was sign-blind). Fits:
  (a) R <= 6 branch (10 pts, all +) with crossing and non-crossing forms;
  (b) |obj| over all 13 with a single-signed form A*G + c*erf(aR)/R.
Same pre-declared bar: rms <= 2e-3 Eh."""
import json
import math
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from f2_stretch import scaled_S, PZPZ  # noqa: E402

d = json.load(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "data", "f2-tail.json")))
pts = [(rec["R"], rec["object"]["pzpz"]) for rec in d["closure_v2"]]
Rall = np.array([p[0] for p in pts])
Yall = np.array([p[1] for p in pts])
kp = 1.1269 ** 2 / 2
Sdiff_all = np.array([float(scaled_S(R, {(0, 0): kp, (1, 0): kp, (0, 1): kp,
                                         (1, 1): kp})[PZPZ]) for R in Rall])


def lin2(g1, g2, y):
    A = np.vstack([g1, g2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef, float(np.sqrt(np.mean((A @ coef - y) ** 2)))


def scan(G, gamfn, Rv, y, sgn=-1.0):
    top = None
    for a in np.arange(0.02, 2.51, 0.01):
        gam = gamfn(Rv, a)
        (A_, c_), rms = lin2(G, sgn * gam, y)
        if top is None or rms < top[0]:
            top = (rms, A_, c_, a)
    return top


erf_over_R = lambda R, a: np.array([math.erf(a * x) for x in R]) / R  # noqa: E731

# ---- (a) positive branch R <= 6
m = Rall <= 6.0
Rb, Yb, Sb = Rall[m], Yall[m], Sdiff_all[m]
print("(a) R<=6 branch, 10 pts:")
for name, G in (("A*S~(0.635) - c*erf/R  ", Sb),
                ("A*S~*(1+.18R)^2 - c*erf/R", Sb * (1 + 0.18 * Rb) ** 2)):
    rms, A_, c_, a = scan(G, erf_over_R, Rb, Yb, sgn=-1.0)
    # where would it cross?
    x = None
    for Rx in np.arange(2.0, 20.01, 0.05):
        Gx = float(scaled_S(Rx, {(0, 0): kp, (1, 0): kp, (0, 1): kp, (1, 1): kp})[PZPZ])
        if "(1+" in name:
            Gx *= (1 + 0.18 * Rx) ** 2
        v = A_ * Gx - c_ * math.erf(a * Rx) / Rx
        if v < 0:
            x = Rx
            break
    print(f"  {name} rms {rms:.2e}  A {A_:+.4f}  c {c_:+.4f}  a {a:.2f}  crossing ~R={x}")

# ---- (b) magnitudes, single-signed
print("(b) |obj| over all 13, single-signed A*G + c*erf/R (A, c same sign):")
Yabs = np.abs(Yall)
for name, G in (("A*|S~(0.635)| + c*erf/R", np.abs(Sdiff_all)),):
    rms, A_, c_, a = scan(G, erf_over_R, Rall, Yabs, sgn=+1.0)
    print(f"  {name} rms {rms:.2e}  A {A_:+.4f}  c {c_:+.4f}  a {a:.2f}")

# ---- (c) how good is the pure-tail claim on the last 4 magnitudes?
tailR, tailY = Rall[-4:], np.abs(Yall[-4:])
for a in (0.2, 0.3, 0.5, 1.0):
    gam = erf_over_R(tailR, a)
    c = float(gam @ tailY / (gam @ gam))
    rms = float(np.sqrt(np.mean((c * gam - tailY) ** 2)))
    print(f"(c) tail-only erf(a={a})/R: c {c:+.4f}  rms {rms:.2e}")
print(f"(c) tail-only c/R: c {float((1/tailR) @ tailY / ((1/tailR) @ (1/tailR))):+.4f}  "
      f"rms {float(np.sqrt(np.mean(((1/tailR) * ((1/tailR) @ tailY / ((1/tailR) @ (1/tailR))) - tailY) ** 2))):.2e}")
