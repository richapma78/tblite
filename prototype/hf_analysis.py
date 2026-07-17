"""hf_analysis.py -- the heteronuclear readings from the HF stretch (offline).

  A. PARAMETER-FREE weight tests (no metric model): Eq. 64's arithmetic mean demands the
     two centers' level responses be EQUAL on every off-diagonal element:
         resp_L2(H,s)/resp_L2(F,s) = 1 on (Hs,Fs);  resp_L2(H,s)/resp_L2(F,p) = 1 on (Hs,Fpz)
     and the mu-law demands the mu responses split 1/2:1/2 x f(q) factors
     (f(+0.116)/f(-0.116) ~ 1.004 from the banked odd-erf samples).
  B. THE FORWARD: H0 = amp * h_bar * Pi_HF * S~sc(mixed metric: k(H)=1.110, kb(H)=0.179,
     k_s(F)=0.730, kb(F)=0.063 -- every number from the H2/F2 fits, NOTHING refit here) with
     amp from named slots (kdiat_sigma(HF) = harmonic(2.909, 2.109) = 2.446), plus the mu
     channel taken DIRECTLY from the measured L7 responses (Euler content, no model).
     Pi_HF = (1 + b_H*R) one-sided with b_H = b_H2 * 2*rcov_H/(rcov_H + rcov_F) -- the
     adapted radii are UNKNOWN (Pyykkoe 0.32/0.64 A assumed, band noted).
  C. THE RESIDUAL = EHT - H0 - mu-channel: the first direct measurement of the POLAR ES
     off-diagonal (+ any q-driven object). Charted against S*q (Mulliken ES1-form leftovers)
     and the L6 (Hubbard) response contents (the ES2 knobs, measured in the collection).
No gate here -- this is a measurement round; the numbers get banked and the ES decode designed
from them.
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

P = params.parse()
eH, eF = P["element"][1], P["element"][9]
kW_s, kW_p = P["globals"]["g1"][0], P["globals"]["g1"][1]
kdiat = 2 * (eH["l8"][0] * eF["l8"][0]) / (eH["l8"][0] + eF["l8"][0])
K_H, KB_H, K_FS, KB_F, B_H2 = 1.110, 0.179, 0.730, 0.063, 0.0315
RCOV_H, RCOV_F = 0.32, 0.64                       # Pyykkoe, Angstrom (adapted set unknown)
B_HF = B_H2 * 2 * RCOV_H / (RCOV_H + RCOV_F)

Bq = basisq.parse()


def s_ham_hf(R):
    """Mixed-metric Ham overlap for HF at separation R (Bohr)."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    cns = adapt.basis_cn([1, 9], xyz)
    sh = []
    for at, (z, k_l, kb) in enumerate((((1), {0: K_H}, KB_H), ((9), {0: K_FS, 1: 1.205}, KB_F))):
        q = kb * math.sqrt(cns[at])
        for l, prims in Bq[z]["shells"]:
            e = np.array([p[0] for p in prims])
            c0 = np.array([p[1] for p in prims])
            c1 = np.array([p[2] for p in prims])
            sh.append({"l": l, "at": at, "exp": e * k_l[l], "coef": (c0 + c1 * q).copy()})
    return overlap.overlap([1, 9], xyz, shells=sh)


def main():
    rows = json.load(open(os.path.join(HERE, "data", "hf-stretch.json")))["rows"]
    print("A -- parameter-free weight tests (Eq. 64 arithmetic mean => ratio = 1):")
    print("   R     L2H/L2F (HsFs)   L2H/L2F (HsFpz)   L7H/L7F (HsFs)")
    for r in rows:
        w1 = r["resp"]["L2_Hs"]["HsFs"] / r["resp"]["L2_Fs"]["HsFs"]
        w2 = r["resp"]["L2_Hs"]["HsFpz"] / r["resp"]["L2_Fp"]["HsFpz"]
        m1 = r["resp"]["L7_Hs"]["HsFs"] / r["resp"]["L7_Fs"]["HsFs"]
        print(f"  {r['R']:5.2f}   {w1:+9.4f}        {w2:+9.4f}         {m1:+9.4f}")

    print("\nB/C -- forward and residual (mixed metric, nothing refit):")
    ampss = kW_s * kdiat
    ampsp = (kW_s + kW_p) / 2 * kdiat
    hss = -(eH["shells"][0][0] + eF["shells"][0][0]) / 2
    hsp = -(eH["shells"][0][0] + eF["shells"][0][1]) / 2
    muH, muFs, muFp = eH["shells"][5][0], eF["shells"][5][0], eF["shells"][5][1]
    print(f"  kdiat_sigma(HF) = {kdiat:.4f}; amp(HsFs) = {ampss:.4f}; amp(HsFpz) = {ampsp:.4f}"
          f"; b_HF = {B_HF:+.4f} (radii band noted)")
    print("   R     EHT(HsFs)  fwd      resid  |  EHT(HsFpz)  fwd      resid  |  S*qH")
    res = []
    for r in rows:
        R = r["R"]
        M = s_ham_hf(R)
        Pi = 1 + B_HF * R
        f1 = ampss * hss * Pi * float(M[0, 1])
        f2 = ampsp * hsp * Pi * float(M[0, 4])
        mu1 = muH * r["resp"]["L7_Hs"]["HsFs"] + muFs * r["resp"]["L7_Fs"]["HsFs"] \
            + muFp * r["resp"]["L7_Fp"]["HsFs"]
        mu2 = muH * r["resp"]["L7_Hs"]["HsFpz"] + muFs * r["resp"]["L7_Fs"]["HsFpz"] \
            + muFp * r["resp"]["L7_Fp"]["HsFpz"]
        r1 = r["EHT"]["HsFs"] - f1 - mu1
        r2 = r["EHT"]["HsFpz"] - f2 - mu2
        res.append((R, r1, r2, r["S"]["HsFs"] * r["q_H"], r["S"]["HsFpz"] * r["q_H"]))
        print(f"  {R:5.2f}  {r['EHT']['HsFs']:+9.5f} {f1 + mu1:+9.5f} {r1:+8.5f} | "
              f"{r['EHT']['HsFpz']:+9.5f} {f2 + mu2:+9.5f} {r2:+8.5f} | {res[-1][3]:+7.4f}")
    print("\nC -- residual correlates (resid / (S*q_H)):")
    for R, r1, r2, sq1, sq2 in res:
        print(f"  {R:5.2f}  HsFs {r1 / sq1 if abs(sq1) > 1e-6 else float('nan'):+8.3f}   "
              f"HsFpz {r2 / sq2 if abs(sq2) > 1e-6 else float('nan'):+8.3f}")
    print("\nL6 (Hubbard/ES2) contents for comparison (slot*resp, HsFs):")
    U_s, U_p = eF["shells"][4][0], eF["shells"][4][1]
    for r in rows:
        c = U_s * r["resp"]["L6_Fs"]["HsFs"] + U_p * r["resp"]["L6_Fp"]["HsFs"]
        print(f"  {r['R']:5.2f}  {c:+9.5f}")


if __name__ == "__main__":
    main()
