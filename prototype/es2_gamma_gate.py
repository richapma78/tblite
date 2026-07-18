"""es2_gamma_gate.py -- regression: the generated ES2 second-order gamma2 must reproduce the
binary's own gamma2 matrix (setespot_ 0x53ab60 arg0, packed) bit-for-bit.

Decoded (NOT fitted) by extracting the real matrix and sweeping the geometry to watch U2(CN):
  gamma2_lAlB = 1 / [R_AB + 0.5(1/U2_lA + 1/U2_lB) * exp(-k2x * R_AB)]        (SI Eq 101)
  U2_lA = T32_lA * ipse_A * (1 + Gamma_A * CN_A)                               (SI Eq 102 form)
  k2x = g2[1] = 0.3300126723 (global). T32=paramfile shell[4], ipse=data/gxtb_ipse.json,
  Gamma=row0[6]. The onsite (R=0) is the harmonic mean of the two U2 -- no separate SRULE.

Reference = the binary's arg0 for HF at 1.733 bohr (CN 0.3062), gdb-extracted via binprobe.
Needs no binary to run.
"""
import math

import paramfile

# binary setespot_ arg0 (packed lower-triangular gamma2) for HF: Hs, Fs, Fp
REF = {"HsHs": 0.5940601242, "FsHs": 0.4100876720, "FsFs": 1.2247713387,
       "FpHs": 0.3726969019, "FpFs": 0.8000639357, "FpFp": 0.5940636595}
CN, SHELLS = 0.3061652499, [("Hs", 0, 1, 0), ("Fs", 1, 9, 0), ("Fp", 1, 9, 1)]
RAT = [[0.0, 1.733], [1.733, 0.0]]


def run():
    P = paramfile.load()
    k2x = P["g2"][1]

    def U2(z, l):
        e = P["elements"][z]
        return e["T32"][l] * e["ipse"] * (1 + e["gamma"] * CN)

    worst = 0.0
    for i in range(3):
        for j in range(i + 1):
            ni, ai, zi, li = SHELLS[i]
            nj, aj, zj, lj = SHELLS[j]
            ua, ub, R = U2(zi, li), U2(zj, lj), RAT[ai][aj]
            g = 1.0 / (R + 0.5 * (1 / ua + 1 / ub) * math.exp(-k2x * R))
            r = REF[ni + nj]
            worst = max(worst, abs(g - r))
            print(f"  {ni}{nj}: gen {g:.10f}  binary {r:.10f}  d {g-r:+.1e}")
    ok = worst < 1e-9
    print(f"\n  worst |generated - binary gamma2| = {worst:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() else 1)
