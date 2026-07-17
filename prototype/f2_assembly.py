"""f2_assembly.py -- the F2 STATIC Fock assembly test: every piece built at the oracle's
converged density, compared element-by-element against the reconstruction.

Pieces:
  H0    forward (file amplitudes x per-shell metric (k_s=0.730, k_p=1.205, kb=0.0632);
        Pi ~ 1; diag = -L2_l) + the L3*CN channel from the MEASURED responses
  ES1   -1/2 S o (v+v), v per shell = mu_l * corr_l(R) with corr extracted from the
        MEASURED mu responses (f2-channels)
  ES2   shell-resolved Mulliken shifts from the restart-P shell charges with the DECODED
        onsite kernel (s-rule x harmonic U) and the scan-level offsite kernel (k2x ~ 0.6/
        labeled HALF-DECODED)
  ACP   analytic matrix
  X     the measured exchange (Euler + the object's closed form on pzpz) -- an R-dependent
        lookup at this stage, NOT a P-functional (labeled)
Residual = F_recon - sum(pieces), reported per element class. This localizes every
assembly error before any SCF loop runs.
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
import constants as K  # noqa: E402
import f2_stretch  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402
import overlap  # noqa: E402

P_ = params.parse()
eF = P_["element"][9]
kW_s, kW_p = P_["globals"]["g1"][0], P_["globals"]["g1"][1]
kdiat_sg, kdiat_pi = eF["l8"][0], eF["l8"][1]
L2s, L2p = eF["shells"][0][0], eF["shells"][0][1]
L3s, L3p = eF["shells"][1][0], eF["shells"][1][1]
MUs, MUp = eF["shells"][5][0], eF["shells"][5][1]
Us, Up = eF["shells"][4][0], eF["shells"][4][1]
K1CN = eF["l1"][8]
KS, KP, KB = 0.730, 1.205, 0.0632
S_F = 0.04193 * 7 + 0.1378          # the s-rule, period 3? F is period 2: a=0.08247 b=0.0920
S_F = 0.08247 * 7 + 0.0920           # = 0.6693

Bq = basisq.parse()
_sh = [(l, np.array([p[0] for p in prims]), np.array([p[1] for p in prims]),
        np.array([p[2] for p in prims])) for l, prims in Bq[9]["shells"]]


def s_ham(R):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    cn = adapt.basis_cn([9, 9], xyz)[0]
    q = KB * math.sqrt(cn)
    kmap = {0: KS, 1: KP}
    sh = [{"l": l, "at": at, "exp": e * kmap[l], "coef": (c0 + c1 * q).copy()}
          for at in (0, 1) for l, e, c0, c1 in _sh]
    return overlap.overlap([9, 9], xyz, shells=sh)


def main():
    ch = json.load(open(os.path.join(HERE, "data", "f2-channels.json")))
    f2rows = {r["R"]: r for r in
              json.load(open(os.path.join(HERE, "data", "f2-stretch.json")))["rows"]}
    obj = json.load(open(os.path.join(HERE, "data", "derived-constants.json")))[
        "hamiltonian_anchors"]["pzpz_object"]["closed_form_MEASURED"]
    for R in (2.668, 3.8):
        rec = fock_recon.fock_ao([("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
        F, S, n = rec["F"], rec["S"], 8
        Sh = s_ham(R)
        # ---- H0
        H0 = np.zeros((n, n))
        lvl = [-L2s, -L2p, -L2p, -L2p] * 2
        for i in range(n):
            H0[i, i] = lvl[i]
        amp = {(0, 0): kW_s * kdiat_sg, (0, 1): (kW_s + kW_p) / 2 * kdiat_sg,
               (1, 1): None}
        for i in range(4):
            for j in range(4, 8):
                li = 0 if i == 0 else 1
                lj = 0 if j == 4 else 1
                sig = (i in (0, 3)) and (j in (4, 7))
                if li == 0 and lj == 0:
                    a = kW_s * kdiat_sg
                elif li != lj:
                    a = (kW_s + kW_p) / 2 * kdiat_sg
                else:
                    a = kW_p * (kdiat_sg if sig else kdiat_pi)
                h = (lvl[i] + lvl[j]) / 2
                H0[i, j] = H0[j, i] = a * h * float(Sh[i, j])
        # L3*CN channel (measured Euler content, element-resolved) -- fold into H0
        crec = ch[str(R)]
        CN_ch = np.zeros((n, n))
        ELS = {"ss": (0, 4), "spz": (0, 7), "pzpz": (3, 7), "pxpx": (1, 5)}
        for el, (i, j) in ELS.items():
            v = (L3s * crec["L3_s"][el] + L3p * crec["L3_p"][el]
                 + K1CN * crec["L1[8]k1cn"][el])
            CN_ch[i, j] = CN_ch[j, i] = v
        # ---- ES1 from measured mu responses (corr per shell from the ll elements)
        corr_s = -(MUs * crec["L7_s"]["ss"] + MUp * crec["L7_p"]["ss"]) / \
            (f2rows[R]["S_ss"] * MUs) if False else None
        # simpler: measured Euler content per element (exact at this geometry)
        ES1 = np.zeros((n, n))
        for el, (i, j) in ELS.items():
            v = MUs * crec["L7_s"][el] + MUp * crec["L7_p"][el]
            ES1[i, j] = ES1[j, i] = v
        # ---- ACP
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        A = f2_stretch.acp_matrix([9, 9], xyz)
        # ---- X measured (Euler) + the object on pzpz
        Xm = np.zeros((n, n))
        for el, (i, j) in ELS.items():
            Xm[i, j] = Xm[j, i] = f2rows[R]["X_el"][el]
        Xm[3, 7] += obj["c"] * math.erf(obj["a"] * R) / R
        Xm[7, 3] = Xm[3, 7]
        # (the object's S~-part folds into the H0-metric side; omitted here -- noted)
        # ---- residual on the four tracked cross elements
        print(f"R = {R}:")
        print("  el      F_recon    H0        CNch      ES1       ACP       Xm        resid")
        for el, (i, j) in ELS.items():
            tot = H0[i, j] + CN_ch[i, j] + ES1[i, j] + A[i, j] + Xm[i, j]
            print(f"  {el:6s} {F[i, j]:+.5f}  {H0[i, j]:+.5f}  {CN_ch[i, j]:+.5f}  "
                  f"{ES1[i, j]:+.5f}  {A[i, j]:+.5f}  {Xm[i, j]:+.5f}  "
                  f"{F[i, j] - tot:+.5f}")
        # diagonals (s and pz on atom 0): assembled = lvl + ES1-diag(?) + ACP + X-diag(?)
        # -- the diagonal channels (mu/L3 responses, X) were measured for CROSS elements
        # only; diagonal assembly needs the atom-gate machinery + onsite ES2: REPORT the
        # raw diagonal residual against (lvl + ACP) as the honest gap.
        for lab, i in (("diag_s", 0), ("diag_pz", 3), ("diag_px", 1)):
            print(f"  {lab:6s} {F[i, i]:+.5f}  vs lvl+ACP {lvl[i] + A[i, i]:+.5f}  "
                  f"gap {F[i, i] - lvl[i] - A[i, i]:+.5f}  (mu/X/ES2-diag not yet assembled)")


if __name__ == "__main__":
    main()
