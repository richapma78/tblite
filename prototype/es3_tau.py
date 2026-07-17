"""tau(R) measured DIRECTLY on H2+ (H-only cation; UKS). For the 2-atom s-only
symmetric case E3_off = q^3 * tau(R) * Gamma_H / 3 exactly (q = +1/2 per atom), so
tau(R) = 3 * clean / (q^3 * Gamma_H) with clean = Gamma_H*dES23/dGamma_H - onsite.
The measured curve is the deliverable; every candidate closed form is then fitted TO
the curve and judged there."""
import math
import re
import sys

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import oracle  # noqa: E402
import params  # noqa: E402

B2A = 1 / 1.8897259886
G_H, U_H = 0.8141717488, 1.0062
K3GS = -0.05708
DELTA = 0.05


def run(R, dG=0.0):
    if dG:
        params.write_perturbed({(1, 1, 6): G_H + dG})
    else:
        params.write_perturbed({})
    r = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R * B2A)], charge=1, uhf=1)
    m = re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)", r["raw"], re.M)
    pops = r.get("pops")
    return float(m.group(1)), pops


curve = []
try:
    for R in (1.2, 1.5, 1.8, 2.1, 2.5, 3.0):
        ep, _ = run(R, +DELTA)
        em, _ = run(R, -DELTA)
        _, pops = run(R, 0.0)
        q = 1.0 - pops[0]["s"] if pops else 0.5
        part = G_H * (ep - em) / (2 * DELTA)
        tau_on = -1.0 / (2 * U_H * U_H)
        onsite = 2 * (1 / 3) * q ** 3 * tau_on * (K3GS * G_H)
        clean = part - onsite
        tau = 3 * clean / (q ** 3 * G_H)
        curve.append((R, q, part, clean, tau))
        print(f"  R={R}: q={q:+.4f}  E3_H-part {part:+.6f}  offsite {clean:+.6f}  "
              f"tau(R) = {tau:+.5f}", flush=True)
finally:
    params.write_perturbed({})

print("\nmeasured tau(R) curve (H-H pair, U=1.0062):")
for R, q, part, clean, tau in curve:
    print(f"  {R:4.1f}  {tau:+.5f}")

# candidate forms fitted TO the curve
import numpy as np  # noqa: E402
Rv = np.array([c[0] for c in curve])
tv = np.array([c[4] for c in curve])
iu = 1.0 / U_H          # invU mean for the H-H pair


def form(k3, k3x, wm):
    w = {"invU": iu, "invU2": iu * iu, "ub": U_H, "ub2": U_H * U_H}[wm]
    return k3 * U_H * Rv * np.exp(-k3x * w * Rv) * (1 - k3x * w * Rv)


best = None
for wm in ("invU", "invU2", "ub", "ub2"):
    for k3 in np.arange(-10, 10.01, 0.1):
        if abs(k3) < 1e-9:
            continue
        for k3x in np.arange(0.02, 3.01, 0.02):
            r = form(k3, k3x, wm) - tv
            r2 = float(r @ r)
            if best is None or r2 < best[0]:
                best = (r2, wm, k3, k3x)
r2, wm, k3, k3x = best
print(f"\nbest form on the curve: w={wm} k3={k3:+.3f} k3x={k3x:.3f}  "
      f"rms {math.sqrt(r2 / len(curve)):.5f}")
for (R, q, part, clean, tau), pred in zip(curve, form(k3, k3x, wm)):
    print(f"  R={R}: tau {tau:+.5f}  form {pred:+.5f}  d {pred - tau:+.5f}")
