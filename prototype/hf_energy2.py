"""hf_energy2.py -- the HF energy re-gate with the three polar wirings.

  ES2+3 offsite: the plain Klopman-Ohno kernel gamma = 1/(R + (1/U_l + 1/U_l')/2)
                 (k2x refit = 0.00 on the banked HCl data; scan-grade, labeled)
  ES1 offsite (Eq 86): -sum drho0_lB * gamma * q_lA (+ mirror); drho0 from REFOCC
  ES1 mu-CN (Eq 84): H-side from the BANKED x-curve (measured 1.5-1.9 Bohr, the 4th-table
                 campaign); F-side calibrated at r_e (IN-SAMPLE, labeled) -- at R=2.7 the
                 CN coupling is dead, so that gate is clean out-of-sample.
GATES: |ES1 d| <= 0.005 and |ES2+3 d| <= 0.01 at R = 2.7 (out-of-sample); r_e reported
with the calibrated factor labeled.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import hf_scf  # noqa: E402

XCURVE = {1.5: 0.3662, 1.7: 0.2164, 1.9: 0.1103}     # banked dES1/dk1H/(mu*q), HF pair
K1H = 0.7749747393                                    # L1[8](H)
DRHO_F = {0: 2.0 - K.REFOCC[9][0], 1: 5.0 - K.REFOCC[9][1]}


def x_h(R):
    xs = sorted(XCURVE)
    if R <= xs[0]:
        return XCURVE[xs[0]]
    if R >= xs[-1]:
        # erf-CN tail: dead by ~2.5 (rc = 0.99); linear-extrapolate to zero at 2.4
        v19 = XCURVE[1.9]
        return max(0.0, v19 * (2.4 - R) / 0.5)
    for a, b in zip(xs, xs[1:]):
        if a <= R <= b:
            t = (R - a) / (b - a)
            return XCURVE[a] * (1 - t) + XCURVE[b] * t


def gamma_ko(R, Ua, Ub):
    return 1.0 / (R + 0.5 * (1.0 / Ua + 1.0 / Ub))


def term(raw, name):
    m = re.search(rf"^\s*{re.escape(name)}\s*:\s*(-?\d+\.\d+)", raw, re.M)
    return float(m.group(1))


def main(cal_F=None):
    results = {}
    for R in (2.7, 1.733):
        rec = fock_recon.fock_ao([("H", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
        P = rec["state"]["P"]
        S = rec["S"]
        raw = rec["state"]["raw"]
        m = np.diag((P / 2.0) @ S)
        qH = 1.0 - 2 * m[0]
        qFs = K.REFOCC[9][0] - 2 * m[1]
        qFp = K.REFOCC[9][1] - 2 * (m[2] + m[3] + m[4])
        qF = qFs + qFp
        muH = hf_scf.MU[("H", 0)]
        muFs, muFp = hf_scf.MU[("F", 0)], hf_scf.MU[("F", 1)]
        uH = hf_scf.UH[("H", 0)]
        uFs, uFp = hf_scf.UH[("F", 0)], hf_scf.UH[("F", 1)]
        fH, fF = 1 + 0.0165 * qH, 1 + 0.0165 * qF
        # ---- ES1: onsite with mu-CN (H measured; F calibrated factor cal_F * x_h shape)
        xh = x_h(R)
        muH_eff = muH * (1 + K1H * xh)
        xf = (cal_F if cal_F is not None else 0.0) * xh
        es1_on = muH_eff * fH * qH + (muFs * qFs + muFp * qFp) * (1 + xf) * fF
        # ---- ES1 offsite (Eq 86)
        es1_off = -(qH * (DRHO_F[0] * gamma_ko(R, uH, uFs)
                          + DRHO_F[1] * gamma_ko(R, uH, uFp)))
        es1 = es1_on + es1_off
        # ---- ES2 onsite + offsite (KO kernel)
        es2 = 0.5 * qH * qH * hf_scf.SRULE["H"] * uH
        for a, ua, qa in ((0, uFs, qFs), (1, uFp, qFp)):
            for b, ub, qb in ((0, uFs, qFs), (1, uFp, qFp)):
                g2 = hf_scf.SRULE["F"] * 2 * ua * ub / (ua + ub)
                es2 += 0.5 * qa * qb * g2
        es2_off = qH * (qFs * gamma_ko(R, uH, uFs) + qFp * gamma_ko(R, uH, uFp))
        es2t = es2 + es2_off
        p1, p23 = term(raw, "ES1 (charge SIE)"), term(raw, "ES2+3")
        results[R] = (es1, p1, es2t, p23, qH, qFs, qFp)
        print(f"R={R}:  q(H {qH:+.4f}; Fs {qFs:+.4f}, Fp {qFp:+.4f})   x_H(R) = {xh:.4f}")
        print(f"  ES1:   ours {es1:+.5f} (on {es1_on:+.5f}, off {es1_off:+.5f})  "
              f"printed {p1:+.5f}  d {es1 - p1:+.5f}")
        print(f"  ES2+3: ours {es2t:+.5f} (on {es2:+.5f}, off {es2_off:+.5f})  "
              f"printed {p23:+.5f}  d {es2t - p23:+.5f}", flush=True)
    es1, p1, *_ = results[2.7]
    es2t = results[2.7][2]
    p23 = results[2.7][3]
    ok = abs(es1 - p1) <= 0.005 and abs(es2t - p23) <= 0.01
    print("OUT-OF-SAMPLE GATE (R=2.7): " + ("PASSED" if ok else "MISSED"))
    return results


if __name__ == "__main__":
    # calibrate the F-side mu-CN factor at r_e: solve cal_F so ES1(1.733) matches
    print("== pass 1 (cal_F = 0) ==")
    r0 = main(0.0)
    es1_0, p1_0 = r0[1.733][0], r0[1.733][1]
    # sensitivity: d(es1)/d(cal_F) = (muFs qFs + muFp qFp) * x_h(1.733) * fF
    qFs, qFp = r0[1.733][5], r0[1.733][6]
    sens = (hf_scf.MU[("F", 0)] * qFs + hf_scf.MU[("F", 1)] * qFp) * 0.199
    cal = (p1_0 - es1_0) / sens
    print(f"\ncalibrated F-side factor = {cal:+.4f} (IN-SAMPLE at r_e, labeled)")
    print("== pass 2 (calibrated) ==")
    main(cal)
