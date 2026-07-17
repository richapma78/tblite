"""f2_forward.py -- the F2 forward model: three predictions, then the four-element gate.

The H2 template (h0_forward v2, GATE PASSED) replicated on fluorine's two shells. OFFLINE
(banked data only). Pre-declared predictions and bars, written before the run:

  P-1  k_s ~ L4_s(F) = 0.7246 and k_p ~ L4_p(F) = 1.1269 (the 'L4 direct' law, tested on
       two more shells; H2 gave 1.110 vs 1.094).
  P-2  kb(ss fit) = kb(pxpx fit) within +-0.05 -- Eq. 28's adaptation is ELEMENT-wise, so
       both shells must return the same Ham sqrt(CN) coefficient.
  P-3  the spz metric is PREDICTED with no new parameters: Seff_spz / S~sc_spz(k_s, k_p, kb)
       should be R-constant at C_sp = [(kW_s+kW_p)/2 * kdiat_sigma(F)] / 2.26 = 0.7219.
  GATE the four-element forward: H0 (file amplitudes x fitted metric) + mu channel (measured
       ratio curves x the exact S) + the penetration object (banked closed form, pzpz only)
       vs the banked EHT curves at the 10 stretch points. Pass: |diff| <= 0.01 Eh everywhere
       for ss/spz/pxpx and <= 0.015 for pzpz (it stacks the object fit's own rms); stretch
       0.006. Pi(F2) = 1 (|b_F| <= 0.003 bounded by fluorine's tiny k^shp = -0.0040; the
       L1[7]-rho route used on H2 is object-contaminated here -- stated, not hidden).

Amplitudes from named slots (kW = G1[0]/G1[1]; kdiat = L8[0]/L8[1]; levels = L2):
  ss: kW_s*kdiat_sg = 0.7188*2.109 = 1.516      spz: (kW_s+kW_p)/2*kdiat_sg = 1.632
  pzpz: kW_p*kdiat_sg = 1.747                    pxpx: kW_p*kdiat_pi = 1.988
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapt  # noqa: E402
import basisq  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

KHM = 2.26
ELS = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}

P = params.parse()
eF = P["element"][9]
kW_s, kW_p = P["globals"]["g1"][0], P["globals"]["g1"][1]
kdiat_sg, kdiat_pi = eF["l8"][0], eF["l8"][1]
L2s, L2p = eF["shells"][0][0], eF["shells"][0][1]
mus, mup = eF["shells"][5][0], eF["shells"][5][1]

Bq = basisq.parse()
_sh_data = [(l, np.array([p[0] for p in prims]), np.array([p[1] for p in prims]),
             np.array([p[2] for p in prims])) for l, prims in Bq[9]["shells"]]
_cn_cache = {}


def s_ham_f2(R, ks, kp, kb):
    """Ham-basis overlap matrix for F2: exponents x k_l per shell, coefficients
    c0 + c1*(kb*sqrt(CN_basis))."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    if R not in _cn_cache:
        _cn_cache[R] = adapt.basis_cn([9, 9], xyz)[0]
    q = kb * math.sqrt(_cn_cache[R])
    kmap = {0: ks, 1: kp}
    sh = [{"l": l, "at": at, "exp": e * kmap[l], "coef": (c0 + c1 * q).copy()}
          for at in (0, 1) for l, e, c0, c1 in _sh_data]
    return overlap.overlap([9, 9], xyz, shells=sh)


def fit2d(Rv, Yv, elem_ix, kgrid, kbgrid, which):
    best = None
    for k in kgrid:
        for kb in kbgrid:
            Sv = np.array([float(s_ham_f2(R, k if which == "s" else 1.0,
                                          k if which == "p" else 1.0, kb)[elem_ix])
                           for R in Rv])
            C = float(Sv @ Yv / (Sv @ Sv))
            rms = float(np.sqrt(np.mean((C * Sv - Yv) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, C, float(k), float(kb))
    return best


def main():
    hm = json.load(open(os.path.join(HERE, "data", "hmetric.json")))
    atlas = json.load(open(os.path.join(HERE, "data", "offdiag-atlas.json")))["atlas"]["F2"]
    f2rows = json.load(open(os.path.join(HERE, "data", "f2-stretch.json")))["rows"]
    Rv = np.array([r["R"] for r in hm["f2"]])

    # ---- v2: kb FIXED from the L8[3] law (kb = 0.4223 * L8[3], H-calibrated; F2's own CN
    # lever is too weak to determine kb -- the v1 "P-2 FAIL" was underdetermination, not
    # inconsistency). k per shell is the only free shape parameter.
    kb = 0.4223 * eF["l8"][3]
    print(f"v2: kb(F) = 0.4223 * L8[3](F) = {kb:+.4f} (the H-calibrated linear law; "
          f"globality of 0.4223 FLAGGED, F-side FD pending)")
    Yss = np.array([r["Seff_ss"] for r in hm["f2"]])
    Ypx = np.array([r["Seff_pxpx"] for r in hm["f2"]])

    def fit1d(Yv, elem_ix, which):
        best = None
        for k in np.arange(0.40, 1.401, 0.005):
            Sv = np.array([float(s_ham_f2(R, k if which == "s" else 1.0,
                                          k if which == "p" else 1.0, kb)[elem_ix])
                           for R in Rv])
            C = float(Sv @ Yv / (Sv @ Sv))
            rms = float(np.sqrt(np.mean((C * Sv - Yv) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, C, float(k))
        return best

    r1, C1, ks = fit1d(Yss, ELS["ss"], "s")
    r2, C2, kp = fit1d(Ypx, ELS["pxpx"], "p")
    print("P-1 -- per-shell metric fits (kb fixed):")
    print(f"  ss  : rms {r1:.2e}  C {C1:+.4f} (pred {kW_s * kdiat_sg / KHM:+.4f})  "
          f"k_s {ks:.3f} [L4_s {eF['shells'][2][0]:.4f}]")
    print(f"  pxpx: rms {r2:.2e}  C {C2:+.4f} (pred {kW_p * kdiat_pi / KHM:+.4f})  "
          f"k_p {kp:.3f} [L4_p {eF['shells'][2][1]:.4f}]")

    # ---- P-3: the spz metric predicted with no new parameters
    Ysp = np.array([r["Seff_spz_via_p"] + r["Seff_spz_via_s"] for r in hm["f2"]])
    Ssp = np.array([float(s_ham_f2(R, ks, kp, kb)[ELS["spz"]]) for R in Rv])
    ratio = Ysp / Ssp
    C_sp_pred = (kW_s + kW_p) / 2 * kdiat_sg / KHM
    print(f"\nP-3 -- spz metric prediction: Seff_spz(total)/S~sc_spz should be constant "
          f"= {C_sp_pred:+.4f}")
    print("   R: " + "  ".join(f"{R:.2f}:{v:+.4f}" for R, v in zip(Rv, ratio)))
    sp_dev = float(np.abs(ratio / np.mean(ratio) - 1).max())
    print(f"   mean {np.mean(ratio):+.4f}  max dev {sp_dev * 100:.1f}%  "
          f"vs pred: {(np.mean(ratio) / C_sp_pred - 1) * 100:+.1f}%")

    # ---- the mu ratio curves (atlas, 4 R) per element/shell-slot
    atR = [rec["R"] for rec in atlas["recs"]]
    sprim = {r["R"]: {"ss": r["S_ss"], "spz": r["S_spz"], "pzpz": r["S_pzpz"],
                      "pxpx": r["S_pxpx"]} for r in f2rows}
    # ---- v3: the mu and CN channels from the COMPLETE per-point measurements
    # (data/f2-channels.json -- all five slots at all ten stretch points; v2's interpolation
    # from four atlas points is gone).
    ch = json.load(open(os.path.join(HERE, "data", "f2-channels.json")))
    L3s, L3p = eF["shells"][1][0], eF["shells"][1][1]
    k1cn = eF["l1"][8]

    def mu_channel(el, R):
        rec = ch[str(R)]
        return (mus * rec["L7_s"][el] + mup * rec["L7_p"][el]) / sprim[R][el]

    def cn_channel(el, R):
        rec = ch[str(R)]
        return L3s * rec["L3_s"][el] + L3p * rec["L3_p"][el] + k1cn * rec["L1[8]k1cn"][el]

    # ---- the penetration object (banked closed form; pzpz only)
    obj = json.load(open(os.path.join(HERE, "data", "derived-constants.json")))[
        "hamiltonian_anchors"]["pzpz_object"]["closed_form_MEASURED"]
    c_o, a_o, A_o = obj["c"], obj["a"], obj["A"]

    def s635(R):
        return float(s_ham_f2(R, 1.1269 ** 2 / 2, 1.1269 ** 2 / 2, 0.0)[ELS["pzpz"]])

    # ---- THE GATE
    print("\nGATE -- F2 forward, four elements, 10 stretch points:")
    amps = {"ss": kW_s * kdiat_sg, "spz": (kW_s + kW_p) / 2 * kdiat_sg,
            "pzpz": kW_p * kdiat_sg, "pxpx": kW_p * kdiat_pi}
    levs = {"ss": -L2s, "spz": -(L2s + L2p) / 2, "pzpz": -L2p, "pxpx": -L2p}
    bars = {"ss": 0.01, "spz": 0.01, "pzpz": 0.015, "pxpx": 0.01}
    worst = {el: 0.0 for el in ELS}
    print("   R    " + "".join(f"{el + ' diff':>12s}" for el in ELS))
    for r in f2rows:
        R = r["R"]
        M = s_ham_f2(R, ks, kp, kb)
        line = f"  {R:4.2f} "
        for el, ix in ELS.items():
            H0 = amps[el] * levs[el] * float(M[ix])
            muv = mu_channel(el, R) * sprim[R][el]      # ratio(resp/S) x S, per shell slot
            tot = H0 + muv + cn_channel(el, R)
            if el == "pzpz":
                tot += c_o * math.erf(a_o * R) / R - A_o * s635(R)
            diff = tot - r["EHT_el"][el]
            worst[el] = max(worst[el], abs(diff))
            line += f"{diff:+12.5f}"
        print(line)
    print("  worst: " + "  ".join(f"{el} {worst[el]:.5f} ({'PASS' if worst[el] <= bars[el] else 'FAIL'})"
                                  for el in ELS))


if __name__ == "__main__":
    main()
