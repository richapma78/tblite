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
  PART 2  the level metric: fit Seff(R) = C * (1 + b*R)^2 * S~(k) per shell pair.
          PRE-DECLARED CHECKS: C*2.26 should equal kW_l * k^diat_L (H2 ss: 0.7188*2.909 =
          2.091; F2 ss: 0.7188*2.109 = 1.516; F2 pi: 0.8284*2.400 = 1.988); k identifies
          the Ham-basis exponent scale for that shell.
  PART 3  the H2 ss FORWARD GATE: H0(file constants + fitted b,k) + mu-channel(metric from
          part 1) vs the banked EHT12 curve. PRE-DECLARED: pass = |forward - measured| <=
          0.01 Eh at every one of the 13 points (stretch: 0.006, the closure noise).
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

    # ---------------- PART 2: the level metric ---------------------------------------
    print("\nPART 2 -- Seff fits: Seff = C*(1+bR)^2*S~(k)")
    kgrid = np.arange(0.30, 1.301, 0.01)
    bgrid = np.arange(-0.10, 0.501, 0.005)
    Rh = np.array([r["R"] for r in hm["h2"]])
    Yh = np.array([r["Seff_ss"] for r in hm["h2"]])
    rms, C, b, k = fit_seff(Rh, Yh, s_scaled_h2, kgrid, bgrid)
    predC = kW_s * kdiat_sH / KHM
    print(f"  H2 ss:  rms {rms:.2e}  C {C:+.4f} (pred {predC:+.4f})  b {b:+.3f}  k {k:.2f}"
          f"   [L4 {eH['shells'][2][0]:.3f}, L4^2/2 {eH['shells'][2][0] ** 2 / 2:.3f}]")
    h2fit = (C, b, k)
    Rf = np.array([r["R"] for r in hm["f2"]])
    for el, key, pred, sel in (
            ("ss", "Seff_ss", kW_s * kdiat_sF / KHM, lambda M: float(M[0, 4])),
            ("pxpx", "Seff_pxpx", kW_p * kdiat_pF / KHM, lambda M: float(M[1, 5]))):
        Yf = np.array([r[key] for r in hm["f2"]])
        def sf(R, k, _sel=sel, _el=el):
            M = s_scaled_f2(R, k if _el == "ss" else 1.0, k if _el != "ss" else 1.0)
            return _sel(M)
        rms, C, b, k = fit_seff(Rf, Yf, sf, kgrid, bgrid)
        l4 = eF["shells"][2][0 if el == "ss" else 1]
        print(f"  F2 {el:4s}: rms {rms:.2e}  C {C:+.4f} (pred {pred:+.4f})  b {b:+.3f}  "
              f"k {k:.2f}   [L4 {l4:.3f}, L4^2/2 {l4 ** 2 / 2:.3f}]")

    # ---------------- PART 2b: b constrained by the L1[7] response -------------------
    # rho_L17 = 2*(kshp_l/Rcov)*R/(1+bR) (the polynomial's own derivative over Pi), so the
    # banked atlas rho curves FIX b independently: 1/(rho/R) is linear in R with slope
    # b/(2a). H2: b = +0.033 (a = 0.280). F's k^shp = -0.0040 is TINY -> Pi(F2) ~ 1, b ~ 0.
    print("\nPART 2b -- Seff refits with b constrained by the L1[7] response:")
    for tag, Rv, Yv, sfun, bfix in (
            ("H2 ss", Rh, Yh, s_scaled_h2, 0.033),
            ("F2 ss", Rf, np.array([r["Seff_ss"] for r in hm["f2"]]),
             lambda R, k: float(s_scaled_f2(R, k, 1.0)[0, 4]), 0.0),
            ("F2 pxpx", Rf, np.array([r["Seff_pxpx"] for r in hm["f2"]]),
             lambda R, k: float(s_scaled_f2(R, 1.0, k)[1, 5]), 0.0)):
        best = None
        for k in kgrid:
            Sv = np.array([sfun(R, k) for R in Rv])
            g = (1 + bfix * Rv) ** 2 * Sv
            Cf = float(g @ Yv / (g @ g))
            rms = float(np.sqrt(np.mean((Cf * g - Yv) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, Cf, k)
        rms, Cf, k = best
        print(f"  {tag:8s} b={bfix:+.3f} fixed:  rms {rms:.2e}  C {Cf:+.4f}  k {k:.2f}")

    # ---------------- PART 3: the H2 ss forward gate ---------------------------------
    print("\nPART 3 -- H2 ss forward gate (pre-declared: |diff| <= 0.01 Eh at all 13 pts):")
    # mu channel: C_mu(R) measured at the 4 atlas points (its mild drift is the mu-CN
    # coupling on the 4th hidden radius table), interpolated in R -- semi-empirical by
    # construction and labeled so. H0 is fully file constants + the part-2 (b, k) shape.
    atR = [rec["R"] for rec in atlas["H2"]["recs"]]
    cmus = [rec["resp"]["L7_s"]["ss"] / s12[rec["R"]] for rec in atlas["H2"]["recs"]]
    C, b, k = h2fit
    L2, mu = eH["shells"][0][0], eH["shells"][5][0]
    worst = 0.0
    print(f"  mu-channel C_mu(R) at {atR} = " + " ".join(f"{c:+.4f}" for c in cmus) +
          "  (interpolated; the drift IS the mu-CN coupling)")
    print(f"  H0 amp = kW_s*kdiat_s(H) = {kW_s * kdiat_sH:+.4f}; level = -L2 (the L3*CN term"
          f" is dropped: bounded <= 0.005 Eh and only at R = 1.4 and below)")
    print("   R     EHT(meas)   H0(fwd)    mu(fwd)    total     diff")
    for r in h2rows:
        R = r["R"]
        H0 = kW_s * kdiat_sH * (-L2) * (1 + b * R) ** 2 * s_scaled_h2(R, k)
        muv = float(np.interp(R, atR, cmus)) * mu * s12[R]
        tot = H0 + muv
        diff = tot - r["EHT12"]
        worst = max(worst, abs(diff))
        print(f"  {R:4.2f}  {r['EHT12']:+9.5f}  {H0:+9.5f}  {muv:+9.5f}  {tot:+9.5f}  "
              f"{diff:+8.5f}")
    print(f"  worst |diff| = {worst:.5f}  -> " +
          ("GATE PASSED" if worst <= 0.01 else "GATE FAILED (report honestly)"))


if __name__ == "__main__":
    main()
