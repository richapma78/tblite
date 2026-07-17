"""h0_forward.py -- the H0 FORWARD MODEL: predict the off-diagonal from the file, then gate.

Eq. 64 (read exactly from the SI PDF, layout mode):

    H_munu = (kW_lA + kW_lB)/2 * (H_lA + H_lB)/2 * Pi_lAlB(R) * S~sc_munu     (A != B)
    H_mumu = H_l (same-atom blocks diagonal)
    H_l    = h_l - k^CN_l * CN                                                 (Eq. 65)
    Pi     = pi_lA * pi_lB,  pi_l = 1 + kshp_A * kshp_l * R/Rcov               (Eq. 66/67)
    S~sc   = diat-frame scaled overlap (sigma/pi/delta x harmonic k^diat)      (Eq. 31/32)
             built on the HAM basis (shell-wise scaled exponents)

All named: kW_s/kW_p = G1[0]/G1[1]; k^diat = L8[0..2]; kshp_A = L1[7]; levels/CN = L2/L3;
mu rides its own ES1 channel. The measured "EHT" curves (F - X - ACP) contain H0 + the mu
channel (+ the pzpz penetration object, fitted separately). This script runs OFFLINE on the
banked datasets -- no oracle calls:

  PART 1  the mu channel's metric: is resp_mu(R) proportional to the PRIMARY overlap
          (ES1 Fock = (1/2)S(v_A+v_B) Mulliken form) or to the level metric? Decided by
          which ratio is R-constant.
  PART 2  the level metric SHAPE-FREE: b from the L1[7] rho (rho = 2aR/(1+bR), banked),
          then M = Seff/Pi matched against S~(k, kb) -- exponents x k, coefficients
          c0 + c1*(kb*sqrt(CN_basis)) -- the Ham basis's own Eq-28 sqrt(CN) adaptation.
          v1's (C, b, k) fits were DEGENERATE; the kb dimension breaks it (v2 result:
          rms 0.01% of range, sharp unique minimum; k = 1.110 ~ L4(H) DIRECT -- the old
          L4^2/2 was an artifact of fitting the level+mu MIXTURE, and the probe-2 FD with
          one delta could not separate the two functional forms; kb = 0.17 vs the density
          basis's 0.2272; C*2.26 = 2.082 vs kW_s*kdiat_s(H) = 2.091 -- 0.4%).
  PART 3  the H2 ss FORWARD GATE v2: H0 = kW_s*kdiat_s(H)*(-L2)*Pi(b)*S~(k, kb) [file
          amplitude, measured shape] + mu-channel = mu*C_mu(R)*S12 with C_mu MEASURED at 7
          R (data/mu-cmu-h2.json; v1 extrapolated below 1.4). PRE-DECLARED: pass =
          |forward - measured| <= 0.01 Eh at all 13 points (stretch: 0.006). The L3*CN
          level term stays out (bounded <= 0.005 Eh, R <= 1.4 only, internal-CN radius
          not decoded for H2).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

KHM = 2.26                      # hmetric's normalization (cancels in the C predictions)


def s_scaled_h2(R, k):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    shells, _ = overlap.build_shells([1, 1], xyz)
    sc = [{"l": s["l"], "at": s["at"], "exp": s["exp"] * k, "coef": s["coef"].copy()}
          for s in shells]
    return float(overlap.overlap([1, 1], xyz, shells=sc)[0, 1])


def s_scaled_f2(R, ks, kp):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    shells, _ = overlap.build_shells([9, 9], xyz)
    kmap = {0: ks, 1: kp}
    sc = [{"l": s["l"], "at": s["at"], "exp": s["exp"] * kmap[s["l"]], "coef": s["coef"].copy()}
          for s in shells]
    return overlap.overlap([9, 9], xyz, shells=sc)


def fit_seff(Rv, Yv, sfun, kgrid, bgrid):
    """Yv ~ C*(1+bR)^2*S~(k): C analytic per (k,b); return best (rms, C, b, k)."""
    best = None
    for k in kgrid:
        Sv = np.array([sfun(R, k) for R in Rv])
        for b in bgrid:
            g = (1 + b * Rv) ** 2 * Sv
            C = float(g @ Yv / (g @ g))
            rms = float(np.sqrt(np.mean((C * g - Yv) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, C, float(b), float(k))
    return best


def main():
    P = params.parse()
    eH, eF = P["element"][1], P["element"][9]
    kW_s, kW_p = P["globals"]["g1"][0], P["globals"]["g1"][1]
    kdiat_sH, kdiat_sF, kdiat_pF = eH["l8"][0], eF["l8"][0], eF["l8"][1]

    hm = json.load(open(os.path.join(HERE, "data", "hmetric.json")))
    atlas = json.load(open(os.path.join(HERE, "data", "offdiag-atlas.json")))["atlas"]
    h2rows = json.load(open(os.path.join(HERE, "data", "h2-eht-probes.json")))["rows"]

    # ---------------- PART 1: the mu channel's metric --------------------------------
    print("PART 1 -- mu metric (constant column = the metric it rides):")
    s12 = {r["R"]: r["S12"] for r in h2rows}
    print("  H2 ss:   R    resp_mu     /S_prim    /Seff")
    seff_h2 = {r["R"]: r["Seff_ss"] for r in hm["h2"]}
    for rec in atlas["H2"]["recs"]:
        R = rec["R"]
        rm = rec["resp"]["L7_s"]["ss"]
        print(f"        {R:4.1f}  {rm:+9.5f}  {rm / s12[R]:+9.5f}  {rm / seff_h2[R]:+9.5f}")
    f2s = json.load(open(os.path.join(HERE, "data", "f2-stretch.json")))["rows"]
    sprim = {r["R"]: r for r in f2s}
    seff_f2 = {r["R"]: r for r in hm["f2"]}
    print("  F2 ss (L7_s):   R    resp_mu     /S_prim    /Seff")
    for rec in atlas["F2"]["recs"]:
        R = rec["R"]
        rm = rec["resp"]["L7_s"]["ss"]
        print(f"        {R:5.2f}  {rm:+9.5f}  {rm / sprim[R]['S_ss']:+9.5f}  "
              f"{rm / seff_f2[R]['Seff_ss']:+9.5f}")
    print("  F2 pxpx (L7_p):  R    resp_mu     /S_prim    /Seff")
    for rec in atlas["F2"]["recs"]:
        R = rec["R"]
        rm = rec["resp"]["L7_p"]["pxpx"]
        print(f"        {R:5.2f}  {rm:+9.5f}  {rm / sprim[R]['S_pxpx']:+9.5f}  "
              f"{rm / seff_f2[R]['Seff_pxpx']:+9.5f}")

    # ---------------- PART 2: the level metric, shape-free (v2) ----------------------
    print("\nPART 2 -- v2: b from the L1[7] rho, then M = Seff/Pi vs S~(k, kb):")
    Ra = np.array([r["R"] for r in atlas["H2"]["recs"]])
    rho = np.array([r["resp"]["L1[7]hbas"]["ss"] / r["resp"]["L2_s"]["ss"]
                    for r in atlas["H2"]["recs"]])
    inv = Ra / rho
    A2 = np.vstack([np.ones_like(Ra), Ra]).T
    (i0, i1), *_ = np.linalg.lstsq(A2, inv, rcond=None)
    a_shp, b_shp = 1 / (2 * i0), i1 / i0
    print(f"  L1[7]-rho: a = {a_shp:.4f}  b = {b_shp:+.4f}  "
          f"max-resid {np.abs(A2 @ np.array([i0, i1]) - inv).max():.3f}")

    import adapt  # noqa: E402
    import basisq  # noqa: E402
    Bq = basisq.parse()
    prims = Bq[1]["shells"][0][1]
    exps0 = np.array([p[0] for p in prims])
    c0v = np.array([p[1] for p in prims])
    c1v = np.array([p[2] for p in prims])

    def s_ham_h2(R, k, kb):
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        cn = adapt.basis_cn([1, 1], xyz)[0]
        coef = c0v + c1v * (kb * np.sqrt(cn))
        sh = [{"l": 0, "at": i, "exp": exps0 * k, "coef": coef.copy()} for i in (0, 1)]
        return float(overlap.overlap([1, 1], xyz, shells=sh)[0, 1])

    Rh = np.array([r["R"] for r in hm["h2"]])
    Yh = np.array([r["Seff_ss"] for r in hm["h2"]])
    M = Yh / (1 + b_shp * Rh) ** 2
    best = None
    for k in np.arange(0.40, 1.301, 0.02):
        for kb in np.arange(-0.40, 0.801, 0.05):
            Sv = np.array([s_ham_h2(R, k, kb) for R in Rh])
            Cx = float(Sv @ M / (Sv @ Sv))
            r_ = float(np.sqrt(np.mean((Cx * Sv - M) ** 2)))
            if best is None or r_ < best[0]:
                best = (r_, Cx, float(k), float(kb))
    rms, C, k, kb = best
    for k2 in np.arange(k - 0.02, k + 0.0201, 0.005):
        for kb2 in np.arange(kb - 0.05, kb + 0.0501, 0.01):
            Sv = np.array([s_ham_h2(R, k2, kb2) for R in Rh])
            Cx = float(Sv @ M / (Sv @ Sv))
            r_ = float(np.sqrt(np.mean((Cx * Sv - M) ** 2)))
            if r_ < rms:
                rms, C, k, kb = r_, Cx, float(k2), float(kb2)
    print(f"  best: rms {rms:.2e} ({rms / np.abs(M).max() * 100:.2f}% of max|M|)  C {C:+.4f}"
          f"  k {k:.3f} [L4(H) {eH['shells'][2][0]:.3f}]  kb {kb:+.2f} [density b 0.2272]")
    print(f"  C*{KHM} = {C * KHM:.4f}  vs kW_s*kdiat_s(H) = {kW_s * kdiat_sH:.4f}  "
          f"({(C * KHM / (kW_s * kdiat_sH) - 1) * 100:+.1f}%)")

    # ---------------- PART 3: the H2 ss forward gate v2 ------------------------------
    print("\nPART 3 -- H2 ss forward gate v2 (pre-declared: |diff| <= 0.01 Eh, all 13 pts):")
    cmu_data = json.load(open(os.path.join(HERE, "data", "mu-cmu-h2.json")))
    cR = sorted(float(x) for x in cmu_data)
    cV = [cmu_data[str(x) if str(x) in cmu_data else repr(x)]["C_mu"] for x in cR]
    L2, mu = eH["shells"][0][0], eH["shells"][5][0]
    worst = 0.0
    print(f"  C_mu measured at {cR}")
    print(f"  H0 amp = kW_s*kdiat_s(H) = {kW_s * kdiat_sH:+.4f} (file), shape Pi(b={b_shp:+.4f})"
          f"*S~(k={k:.3f}, kb={kb:+.2f}); level = -L2; L3*CN dropped (<= 0.005, R <= 1.4)")
    print("   R     EHT(meas)   H0(fwd)    mu(fwd)    total     diff")
    for r in h2rows:
        R = r["R"]
        H0 = kW_s * kdiat_sH * (-L2) * (1 + b_shp * R) ** 2 * s_ham_h2(R, k, kb)
        muv = float(np.interp(R, cR, cV)) * mu * s12[R]
        tot = H0 + muv
        diff = tot - r["EHT12"]
        worst = max(worst, abs(diff))
        print(f"  {R:4.2f}  {r['EHT12']:+9.5f}  {H0:+9.5f}  {muv:+9.5f}  {tot:+9.5f}  "
              f"{diff:+8.5f}")
    print(f"  worst |diff| = {worst:.5f}  -> " +
          ("GATE PASSED" + (" (stretch 0.006 too)" if worst <= 0.006 else "")
           if worst <= 0.01 else "GATE FAILED (report honestly)"))


if __name__ == "__main__":
    main()
