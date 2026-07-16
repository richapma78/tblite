"""offsite_hcl.py -- decode the OFFSITE second-order kernel gamma(2)(R) from an HCl stretch.

SI Eq. 101: gamma_lAlB = 1 / (R + 1/2*(1/U_lA + 1/U_lB) * exp(-k2x * R)), with U = U0*(1 + kU*CN)
(kU element-wise, CN = the internal Eq.47 CN already decoded for repulsion). The SAME kernel
carries the offsite FIRST order (SI Eq. 86): E1_off = -sum drho0_lB * gamma * q_lA with
drho0 = aufbau - refocc (both known!). Two printed observables per geometry -> two independent
probes of one kernel.

Per stretch point (populations read pointwise -- the v2 lesson):
  offsite_ES23(R) = ES23_printed - onsite_model(H) - onsite_model(Cl)
  offsite_ES1(R)  = ES1_printed  - onsite_ES1_model
with onsite models from the charted pieces (s-rule, degree-2 third order, f-switching linearized;
E4 negligible at |q| ~ 0.1). Then fit k2x (+ optionally kU) against BOTH curves jointly.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

BOHR = K.BOHR
S_H, S_CL = 0.4726, 0.4312          # measured s values (derived-constants.json)
C1_CL = 0.0075                       # ~S-row value; E3 at |q|~0.1 is ~1e-5 -- noise here
F_SLOPE = 0.0165                     # f(q) ~ 1 + F_SLOPE*q near 0 (measured samples)

AUFBAU = {1: {0: 1.0}, 17: {0: 2.0, 1: 5.0, 2: 0.0}}


def onsite_es23(q, U, s, c1):
    E2 = s * sum(0.5 * q[a] * q[b] * (2 * U[a] * U[b] / (U[a] + U[b]))
                 for a in range(len(q)) for b in range(len(q)))
    qa = sum(q)
    E3 = c1 * qa * sum(q[a] * q[b] * 0.25 * (U[a] + U[b]) ** 2
                       for a in range(len(q)) for b in range(len(q)))
    return E2 + E3


def main():
    P = params.parse()
    eH, eCl = P["element"][1], P["element"][17]
    U_H = [eH["shells"][4][0]]
    U_CL = eCl["shells"][4][:3]
    muH = [eH["shells"][5][0]]
    muCl = eCl["shells"][5][:3]
    refH, refCl = K.REFOCC[1], K.REFOCC[17]
    drho_H = [AUFBAU[1][0] - refH[0]]
    drho_CL = [AUFBAU[17][l] - refCl[l] for l in range(3)]

    RS = [2.0, 2.2, 2.409, 2.7, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
    rows = []
    for R in RS:
        r = oracle.run([("Cl", 0, 0, 0), ("H", 0, 0, R / BOHR)])
        pCl, pH = r["pops"][0], r["pops"][1]
        qCl = [refCl[0] - pCl["s"], refCl[1] - pCl["p"], refCl[2] - pCl["d"]]
        qH = [refH[0] - pH["s"]]
        es1, es23 = r["terms"]["ES1 (charge SIE)"], r["terms"]["ES2+3"]
        on23 = onsite_es23(qH, U_H, S_H, 0.0) + onsite_es23(qCl, U_CL, S_CL, C1_CL)
        qA_H, qA_CL = sum(qH), sum(qCl)
        on1 = (sum(m * (1 + F_SLOPE * qA_H) * q for m, q in zip(muH, qH))
               + sum(m * (1 + F_SLOPE * qA_CL) * q for m, q in zip(muCl, qCl)))
        rows.append({"R": R, "qH": qH, "qCl": qCl,
                     "off23": es23 - on23, "off1": es1 - on1})
        print(f"  R={R:4.2f}  qH {qH[0]:+.4f}  qCl_s {qCl[0]:+.4f}  "
              f"off23 {es23 - on23:+9.6f}  off1 {es1 - on1:+9.6f}", flush=True)

    # joint fit of k2x (U at U0: CN-dependence absorbed later; start simple)
    def gamma(R, Ua, Ub, k2x):
        return 1.0 / (R + 0.5 * (1.0 / Ua + 1.0 / Ub) * math.exp(-k2x * R))

    def model(row, k2x):
        R = row["R"]
        o23 = sum(row["qH"][a] * row["qCl"][b] * gamma(R, U_H[a], U_CL[b], k2x)
                  for a in range(1) for b in range(3))
        o1 = -(sum(drho_CL[b] * gamma(R, U_H[0], U_CL[b], k2x) * row["qH"][0] for b in range(3))
               + sum(drho_H[0] * gamma(R, U_H[0], U_CL[b], k2x) * row["qCl"][b] for b in range(3)))
        return o23, o1

    print("\n  k2x scan (joint rms over off23+off1):")
    best = None
    for k2x in np.arange(0.0, 3.01, 0.05):
        sq = 0.0
        for row in rows:
            m23, m1 = model(row, k2x)
            sq += (m23 - row["off23"]) ** 2 + (m1 - row["off1"]) ** 2
        rms = math.sqrt(sq / (2 * len(rows)))
        if best is None or rms < best[0]:
            best = (rms, k2x)
    rms, k2x = best
    print(f"  best k2x = {k2x:.2f}   rms = {rms:.2e}")
    for row in rows:
        m23, m1 = model(row, k2x)
        print(f"    R={row['R']:4.2f} off23 {row['off23']:+9.6f}/{m23:+9.6f}  "
              f"off1 {row['off1']:+9.6f}/{m1:+9.6f}")
    json.dump(rows, open(os.path.join(HERE, "data", "offsite-hcl.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
