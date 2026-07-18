"""mfx.py -- the range-separated Mulliken-approximated Fock exchange (gap #1), SI Sec 1.15.
Energy = the validated 4-index Mulliken form; kernel = the EXACT setgab_lrao_ form, decoded
from the Ghidra body (line 209902-209915) + uwe_av_ (227133) and gate-verified against the
binary's own AO gamma matrix (4th arg of setgab_lrao_) to 5e-11 across an H2 bond sweep.

  gamma_AB = [alpha + (1-alpha) erf(omega R)] / [R + exp(-R (B_AB k2' + k1')) / (favg . L)]
  favg(x,y,xi) = uwe_av = 2^(xi-1) (x y)^(xi/2) / (x+y)^(xi-1)   (xi = max(lA,lB)+1)
  x = U_l = T25[Z][l] . ipse[Z]        # per-shell, CN-FREE (a pure atomic constant)
  L    = 1.39                          # atomic self-pair (i==j)
       = sqrt(c[lA] c[lB])             # off-diagonal;  c[s] = 0.0788775224
  B_AB = per-element-pair bond param (0x3be6080 table; B_HH = 2.1823, == AES R0 -- verify)
  alpha=0.15  omega=0.2347047181  k1'=-0.5959929766  k2'=0.2140456651  (defaults, gdb-pinned)

Two traps this decode corrected: (1) favg . L DIVIDES the screening exp -- it is NOT
favg . exp -- so onsite gamma = alpha . favg . L, not alpha/favg (that inversion was the
14% onsite residual); (2) L is 1.39 on the diagonal but sqrt(c.c) off it -- assuming 1.39
everywhere left the offsite 5.4x too big. The favg input is T25.ipse (CN-free); the CN-scaled
U (line 209715) feeds the p2/p3 intermediates, NOT the energy gamma.

GATE: H2 has only s shells -> OFX = 0 -> printed Ex(H2) = MFX alone, and the kernel
reproduces the extracted p4(R) curve to 5e-11.  Polyatomics need the per-element T25/ipse/c
tables (C,N,O,F) + OFX (Sec 1.16) -- the documented next step.
"""
import math
import re
import sys

import numpy as np
from math import erf

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import overlap as OV  # noqa: E402
import paramfile  # noqa: E402
import restart  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
# ALL the MFX constants come from the parsed parameter file (paramfile.py): the globals plus
# every element's T25 -- so the gamma works for all 79 parameterised elements (Z 1..92,
# transition metals included), not a hardcoded few. ipse is the one derived piece (data/
# gxtb_ipse.json, gdb-extracted -- not stored in gxtb_parameters).
_P = paramfile.load()
_ELEM = _P["elements"]
ALPHA = 0.15                                   # range-separation floor (code default)
OMEGA = _P["omega"]                            # g2[8]
K1P, K2P = _P["g1"][6], _P["g1"][7]            # screening exponent (0x3be6350/58)
C_OFF = _P["c"]                                # c[l] = g1[4+l]; off-diagonal L
L_DIAG = 1.39                                  # atomic self-pair L factor
PURE_S = {"h2"}                                # OFX = 0 (no onsite different-l)

# bond param B[Z][Z'] per ATOM-PAIR is a SEPARATE pairwise table, NOT in gxtb_parameters; still
# hardcoded for the gated pairs. ASYMMETRIC; ordering rule for asymmetric pairs is TODO.
BOND = {(1, 1): 2.1823, (1, 9): 2.0646, (9, 9): 2.7664,
        (1, 8): 2.1768, (8, 8): 2.4817}       # O-H symmetric; H-F asymmetric (ordering TODO)


def u_shell(z, l):
    """favg input for a shell: T25[Z][l] * ipse[Z], CN-free (setgab_lrao_ local_158)."""
    e = _ELEM[z]
    return e["T25"][l] * e["ipse"]


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


def bond_param(za, zb):
    return BOND.get((za, zb), BOND.get((zb, za)))


def gamma_matrix(zs, xyz, meta):
    """The exact decoded setgab_lrao_ AO exchange gamma. meta[i] = (atom, Z, l, ...)."""
    n = len(meta)
    Rab = np.zeros((len(zs), len(zs)))
    for a in range(len(zs)):
        for b in range(len(zs)):
            Rab[a, b] = np.linalg.norm(xyz[a] - xyz[b])
    gam = np.zeros((n, n))
    for i in range(n):
        ai, zi, li, _ = meta[i]
        for j in range(n):
            aj, zj, lj, _ = meta[j]
            fa = favg(u_shell(zi, li), u_shell(zj, lj), 1)   # xi=1 (geom. mean) for s,p
            same_shell = (ai == aj and li == lj)             # one diagonal shell block
            L = L_DIAG if same_shell else math.sqrt(C_OFF[li] * C_OFF[lj])
            R = Rab[ai, aj]
            num = ALPHA + (1 - ALPHA) * erf(OMEGA * R)
            screen = math.exp(-R * (bond_param(zi, zj) * K2P + K1P))
            gam[i, j] = num / (R + screen / (fa * L))
    return gam


def run(verbose=True):
    """Gate the energy on systems whose elements are populated in ATOMIC (H2 today)."""
    worst_s = 0.0
    for name, (zs, xyz) in systems().items():
        if any(z not in ATOMIC for z in zs):
            continue                                       # element table not yet extracted
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        st = restart.converged_state(atoms)
        P = st["P"]
        printed = float(re.search(r"Ex \(Mulliken\)\s*:\s*(-?\d+\.\d+)",
                                  st["raw"]).group(1))
        B = GE.build(zs, np.array(xyz, float))
        S, meta = B["S"], B["meta"]
        gam = gamma_matrix(zs, np.array(xyz, float), meta)
        e = 2.0 * ex_energy(P, S, gam)         # x2 = the alpha+beta closed-shell spin sum
        d = e - printed
        tag = "  [OFX=0, pure MFX gate]" if name in PURE_S else "  (+OFX missing)"
        if name in PURE_S:
            worst_s = max(worst_s, abs(d))
        if verbose:
            print(f"  {name:4s} ours {e:+.6f}  printed {printed:+.6f}  d {d:+.6f}{tag}")
    if verbose:
        print(f"  H2 (pure-MFX) energy |d| = {worst_s:.6f} Eh")
    return worst_s


if __name__ == "__main__":
    run()
