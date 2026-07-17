"""mfx.py -- the range-separated Mulliken-approximated Fock exchange (gap #1), SI Sec 1.15.
Energy = the validated 4-index Mulliken form; kernel = Eq 149 with gdb-extracted params.

  gamma_lAlB = [alpha + (1-alpha) erf(omega R)] / [R + favg(U_lA,U_lB,xi) exp(-(k1+k2 R)R)]
  favg(X,Y,xi) = 2^(xi-1) (XY)^(xi/2) / (X+Y)^(xi-1)   (Eq 150; xi=1 valence, 2 polar)
  alpha=0.15, omega=0.2347047181 (defaults, gdb-confirmed at setgab_lrao_);
  U^MFX = the standard ES Hubbard U (gp3_gam2 = elem U), CN-scaled (line 209715).

GATE strategy: H2 has only s shells -> OFX = 0 -> printed Ex(H2) = MFX alone, so H2 gates
the KERNEL directly. Polyatomics with p shells also carry OFX (Sec 1.16), added later.
Sweeps the uncertain pieces (k1/k2 roles, favg xi, CN-scaling) against the H2 gate.
"""
import math
import re
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import overlap as OV  # noqa: E402
import restart  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
ALPHA, OMEGA = 0.15, 0.2347047181
K1, K2 = 0.0788775224, 1.7995847408          # G1[4], G1[5] (candidate screening)
PURE_S = {"h2"}                               # OFX = 0 (no onsite different-l)


def systems():
    ang = 104.5 * math.pi / 180
    r_oh = 0.9572 * BOHR
    W = [[0.0, 0.0, 0.0],
         [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
         [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
    a4 = 1.087 * BOHR / math.sqrt(3)
    r_nh = 1.012 * BOHR
    st, ct = 0.9262, -0.3770
    return {
        "h2": ([1, 1], np.array([[0, 0, 0], [0, 0, 1.4]], float)),
        "f2": ([9, 9], np.array([[0, 0, 0], [0, 0, 2.668]], float)),
        "hf": ([1, 9], np.array([[0, 0, 0], [0, 0, 1.733]], float)),
        "h2o": ([8, 1, 1], np.array(W)),
        "ch4": ([6, 1, 1, 1, 1], np.array(
            [[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4], [-a4, a4, -a4], [-a4, -a4, a4]],
            float)),
        "nh3": ([7, 1, 1, 1], np.array(
            [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                            r_nh * st * math.sin(2 * math.pi * k / 3),
                            r_nh * ct] for k in range(3)], float)),
    }


def favg(X, Y, xi):
    return 2.0 ** (xi - 1) * (X * Y) ** (xi / 2.0) / (X + Y) ** (xi - 1)


def ex_energy(P, S, gam):
    """the validated 4-index Mulliken exchange energy (scf_h2 conventions)."""
    Ph = P / 2.0
    n = P.shape[0]
    E = 0.0
    for mu in range(n):
        for nu in range(n):
            for lam in range(n):
                for kap in range(n):
                    E -= (Ph[mu, nu] * S[mu, lam] * Ph[lam, kap] * S[kap, nu]
                          * (gam[mu, kap] + gam[mu, nu] + gam[lam, kap]
                             + gam[lam, nu]) / 16.0)
    return 2.0 * E


def gamma_matrix(zs, xyz, meta, U, xi_val=1, k1=K1, k2=K2, ucn=None):
    n = len(meta)
    Rab = np.zeros((len(zs), len(zs)))
    for a in range(len(zs)):
        for b in range(len(zs)):
            Rab[a, b] = np.linalg.norm(xyz[a] - xyz[b])
    gam = np.zeros((n, n))
    for i in range(n):
        ai, zi, li, _ = meta[i]
        Ui = (ucn[i] if ucn is not None else U[i])
        for j in range(n):
            aj, zj, lj, _ = meta[j]
            Uj = (ucn[j] if ucn is not None else U[j])
            xi = xi_val
            fa = favg(Ui, Uj, xi)
            if ai == aj:                                   # onsite (R=0)
                gam[i, j] = ALPHA / fa
            else:
                R = Rab[ai, aj]
                num = ALPHA + (1 - ALPHA) * erf(OMEGA * R)
                den = R + fa * math.exp(-(k1 + k2 * R) * R)
                gam[i, j] = num / den
    return gam


def run(k1=K1, k2=K2, xi_val=1, verbose=True):
    worst_s = 0.0
    for name, (zs, xyz) in systems().items():
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        st = restart.converged_state(atoms)
        P = st["P"]
        printed = float(re.search(r"Ex \(Mulliken\)\s*:\s*(-?\d+\.\d+)",
                                  st["raw"]).group(1))
        B = GE.build(zs, np.array(xyz, float))
        S, meta = B["S"], B["meta"]
        # U^MFX per AO = the standard ES Hubbard U (gam2 = elem U)
        U = np.array([B["E"][meta[i][1]]["U"][meta[i][2]] for i in range(len(meta))])
        gam = gamma_matrix(zs, np.array(xyz, float), meta, U, xi_val=xi_val,
                           k1=k1, k2=k2)
        e = ex_energy(P, S, gam)
        d = e - printed
        tag = "  [OFX=0, pure MFX gate]" if name in PURE_S else "  (+OFX missing)"
        if name in PURE_S:
            worst_s = max(worst_s, abs(d))
        if verbose:
            print(f"  {name:4s} ours {e:+.6f}  printed {printed:+.6f}  d {d:+.6f}{tag}")
    if verbose:
        print(f"  H2 (pure-MFX) |d| = {worst_s:.6f}  "
              f"[k1={k1:.4f} k2={k2:.4f} xi={xi_val}]")
    return worst_s


if __name__ == "__main__":
    run()
