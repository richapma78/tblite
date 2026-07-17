"""energy_gate2.py -- the F2 and HF energy gates, per-term resolved.

The printed decomposition allows term-by-term gating:
  ES1   ours = sum_A sum_l mu_l * f(q_A) * q_l   (f ~ 1 + 0.0165q; decoded)
  ES2+3 ours = onsite shell-resolved E2 (decoded kernel) [+ E3 ~ c1*q_A*... tiny, bounded;
        offsite E2 with the HALF-DECODED kernel omitted -> expect a ~0.005 gap, labeled]
  EHT+ACP ours = Tr(H0_fwd P) + Tr(ACP P)   gated against
        printed(electronic - ES total - Ex - Espinpol)
Ex / AES / dispersion: pass-through (the F2/HF exchange kernel matrices await their own
extraction campaign; AES undecoded). GATES (pre-declared): |ES1| d <= 0.005;
|ES2+3| d <= 0.01; |EHT+ACP| d <= 0.03 at two R per molecule.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import f2_scf  # noqa: E402
import f2_stretch  # noqa: E402
import fock_recon  # noqa: E402
import hf_scf  # noqa: E402

F2U = [f2_scf.U[0], f2_scf.U[1]]
SF = f2_scf.S_F


def terms_of(raw):
    out = {}
    for t in ("ES1 (charge SIE)", "ES2+3", "ES multipole", "ES total", "Ex (Mulliken)",
              "Espinpol", "electronic"):
        m = re.search(rf"^\s*{re.escape(t)}\s*:\s*(-?\d+\.\d+)", raw, re.M)
        out[t] = float(m.group(1)) if m else 0.0
    return out


def gate_f2(R):
    rec = fock_recon.fock_ao([("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
    P = rec["state"]["P"]
    S = rec["S"]
    t = terms_of(rec["state"]["raw"])
    m = np.diag((P / 2.0) @ S)
    q = [K.REFOCC[9][0] - 2 * m[0], K.REFOCC[9][1] - 2 * (m[1] + m[2] + m[3])]
    # ES1 (q_A = 0 -> f = 1), both atoms
    es1 = 2 * (f2_scf.MU[0] * q[0] + f2_scf.MU[1] * q[1])
    # ES2 onsite, both atoms
    es2 = 0.0
    for a in range(2):
        for b in range(2):
            g2 = SF * 2 * F2U[a] * F2U[b] / (F2U[a] + F2U[b])
            es2 += 0.5 * q[a] * q[b] * g2
    es2 *= 2
    f2_scf.INCLUDE_SHORT = False
    S_, H0, A, _OBJ = f2_scf.build_static(R)
    f2_scf.INCLUDE_SHORT = True
    e_eht = float(np.sum(H0 * P))
    e_acp = float(np.sum(A * P))
    ref_ehtacp = t["electronic"] - t["ES total"] - t["Ex (Mulliken)"] - t["Espinpol"]
    print(f"F2 R={R}:")
    print(f"  ES1:   ours {es1:+.5f}  printed {t['ES1 (charge SIE)']:+.5f}  "
          f"d {es1 - t['ES1 (charge SIE)']:+.5f}")
    print(f"  ES2+3: ours {es2:+.5f}  printed {t['ES2+3']:+.5f}  "
          f"d {es2 - t['ES2+3']:+.5f}")
    print(f"  EHT+ACP: ours {e_eht + e_acp:+.5f}  ref(elec-ES-Ex-spin) {ref_ehtacp:+.5f}"
          f"  d {e_eht + e_acp - ref_ehtacp:+.5f}", flush=True)
    return (abs(es1 - t["ES1 (charge SIE)"]) <= 0.005 and
            abs(es2 - t["ES2+3"]) <= 0.01 and
            abs(e_eht + e_acp - ref_ehtacp) <= 0.03)


def gate_hf(R):
    rec = fock_recon.fock_ao([("H", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
    P = rec["state"]["P"]
    S = rec["S"]
    t = terms_of(rec["state"]["raw"])
    m = np.diag((P / 2.0) @ S)
    qH = 1.0 - 2 * m[0]
    qFs = K.REFOCC[9][0] - 2 * m[1]
    qFp = K.REFOCC[9][1] - 2 * (m[2] + m[3] + m[4])
    qF = qFs + qFp
    fH, fF = 1 + 0.0165 * qH, 1 + 0.0165 * qF
    es1 = hf_scf.MU[("H", 0)] * fH * qH + \
        (hf_scf.MU[("F", 0)] * qFs + hf_scf.MU[("F", 1)] * qFp) * fF
    es2 = 0.0
    uH = hf_scf.UH[("H", 0)]
    es2 += 0.5 * qH * qH * hf_scf.SRULE["H"] * uH
    uF = [hf_scf.UH[("F", 0)], hf_scf.UH[("F", 1)]]
    qFl = [qFs, qFp]
    for a in range(2):
        for b in range(2):
            g2 = hf_scf.SRULE["F"] * 2 * uF[a] * uF[b] / (uF[a] + uF[b])
            es2 += 0.5 * qFl[a] * qFl[b] * g2
    hf_scf.INCLUDE_RES = False
    S_, H0, A, _R = hf_scf.build_static(R)
    hf_scf.INCLUDE_RES = True
    e_eht = float(np.sum(H0 * P))
    e_acp = float(np.sum(A * P))
    ref = t["electronic"] - t["ES total"] - t["Ex (Mulliken)"] - t["Espinpol"]
    print(f"HF R={R}:")
    print(f"  q: H {qH:+.4f}  F {qF:+.4f}")
    print(f"  ES1:   ours {es1:+.5f}  printed {t['ES1 (charge SIE)']:+.5f}  "
          f"d {es1 - t['ES1 (charge SIE)']:+.5f}")
    print(f"  ES2+3: ours {es2:+.5f}  printed {t['ES2+3']:+.5f}  "
          f"d {es2 - t['ES2+3']:+.5f}   (offsite kernel omitted, labeled)")
    print(f"  EHT+ACP: ours {e_eht + e_acp:+.5f}  ref {ref:+.5f}  "
          f"d {e_eht + e_acp - ref:+.5f}", flush=True)
    return (abs(es1 - t["ES1 (charge SIE)"]) <= 0.005 and
            abs(es2 - t["ES2+3"]) <= 0.01 and
            abs(e_eht + e_acp - ref) <= 0.03)


def main():
    ok = True
    for R in (2.668, 3.8):
        ok &= gate_f2(R)
    for R in (1.733, 2.7):
        ok &= gate_hf(R)
    print("F2/HF ENERGY GATES " + ("PASSED" if ok else "FAILED (per-term above)"))


if __name__ == "__main__":
    main()
