"""multipole.py -- dipole and quadrupole moment integrals over the adapted q-vSZP basis,
the first layer of the AES (anisotropic electrostatics) port. Extends overlap.py's gated
Obara-Saika engine: for a primitive centred at A, x*g_i = g_{i+1} + A_x*g_i, so the
moment 1-D factors are

  T0 = s[a,b]                                   (overlap)
  T1 = s[a+1,b] + A*s[a,b]                       (<i|x|j>, about origin)
  T2 = s[a+2,b] + 2A*s[a+1,b] + A^2*s[a,b]       (<i|x^2|j>)

built from the SAME _os_1d recursion (evaluated to la+2), the SAME PySCF primitive
normalization and spherical transform as the gated overlap -- so S comes out bit-identical
to overlap.overlap() and D/Q ride the same machinery.

GATE (declared): our D (3) and Q (6 unique) equal PySCF int1e_r / int1e_rr for the same
adapted basis, max |d| <= 1e-10, on the neutral gate systems. Origin = (0,0,0), matching
PySCF's default gauge; CAMM shifts to atomic centres later.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as _K  # noqa: E402
import overlap as OV  # noqa: E402

_ORACLE_D_PERM = _K.D_PERM


def _os_1d(la, lb, pa, pb, mu):
    """Obara-Saika 1-D overlap table s[i,j], i in 0..la, j in 0..lb (unnormalized)."""
    s = np.zeros((la + 1, lb + 1))
    s[0, 0] = 1.0
    inv2mu = 0.5 / mu
    for i in range(1, la + 1):
        s[i, 0] = pa * s[i - 1, 0] + (i - 1) * inv2mu * (s[i - 2, 0] if i > 1 else 0.0)
    for j in range(1, lb + 1):
        for i in range(la + 1):
            s[i, j] = (pb * s[i, j - 1] + i * inv2mu * (s[i - 1, j - 1] if i else 0.0)
                       + (j - 1) * inv2mu * (s[i, j - 2] if j > 1 else 0.0))
    return s


# the 6 unique quadrupole components in PySCF int1e_rr's row order (xx,xy,xz,yy,yz,zz)
_Q6 = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def cart_moments(a, b, ra, rb, la, lb):
    """S, D(3), Q(6) blocks over CARTESIAN components for one primitive pair (unnormalized
    prims), moments about the ORIGIN."""
    mu = a + b
    p = (a * ra + b * rb) / mu
    ab2 = float(np.dot(ra - rb, ra - rb))
    pref = math.exp(-a * b / mu * ab2) * (math.pi / mu) ** 1.5
    # 1-D tables to (la+2) so s[a+1], s[a+2] are valid
    tab = [_os_1d(la + 2, lb, p[k] - ra[k], p[k] - rb[k], mu) for k in range(3)]
    ca, cb = OV.cart_components(la), OV.cart_components(lb)
    nca, ncb = len(ca), len(cb)
    S = np.empty((nca, ncb))
    D = np.empty((3, nca, ncb))
    Q = np.empty((6, nca, ncb))
    for i, ai in enumerate(ca):
        for j, bj in enumerate(cb):
            t0 = [tab[k][ai[k], bj[k]] for k in range(3)]
            t1 = [tab[k][ai[k] + 1, bj[k]] + ra[k] * tab[k][ai[k], bj[k]]
                  for k in range(3)]
            t2 = [tab[k][ai[k] + 2, bj[k]] + 2 * ra[k] * tab[k][ai[k] + 1, bj[k]]
                  + ra[k] * ra[k] * tab[k][ai[k], bj[k]] for k in range(3)]
            S[i, j] = pref * t0[0] * t0[1] * t0[2]
            for d in range(3):
                f = [t0[0], t0[1], t0[2]]
                f[d] = t1[d]
                D[d, i, j] = pref * f[0] * f[1] * f[2]
            for qi, (dA, dB) in enumerate(_Q6):
                f = [t0[0], t0[1], t0[2]]
                if dA == dB:
                    f[dA] = t2[dA]
                else:
                    f[dA] = t1[dA]
                    f[dB] = t1[dB]
                Q[qi, i, j] = pref * f[0] * f[1] * f[2]
    return S, D, Q


def moment_matrices(zs, xyz_bohr, charge=0, shells=None, ao_order="oracle"):
    """Spherical-AO S, D(3,nao,nao), Q(6,nao,nao) for the adapted basis, moments about
    the origin. Same normalization + spherical transform as overlap.overlap()."""
    xyz = np.asarray(xyz_bohr, float)
    if shells is None:
        shells, _ = OV.build_shells(zs, xyz, charge=charge)
    mats = {l: OV._c2s(l) for l in sorted({sh["l"] for sh in shells})}
    if ao_order == "oracle":
        if any(sh["l"] >= 3 for sh in shells):
            raise NotImplementedError("oracle f-shell ordering not measured yet")
        if 2 in mats:
            T = np.zeros((5, 5))
            for o, pp in enumerate(_ORACLE_D_PERM):
                T[pp, o] = 1.0
            mats[2] = T @ mats[2]
    for sh in shells:
        c = sh["coef"] * np.array([OV.prim_norm(a, sh["l"]) for a in sh["exp"]])
        l, e = sh["l"], sh["exp"]
        dfac = 1.0
        for k in range(2 * l - 1, 0, -2):
            dfac *= k
        ee = e[:, None] + e[None, :]
        s_rad = (np.pi ** 1.5) * dfac / (2.0 ** l) / ee ** (l + 1.5)
        sh["cn"] = c / math.sqrt(float(c @ s_rad @ c))
    nao = sum(2 * sh["l"] + 1 for sh in shells)
    S = np.zeros((nao, nao))
    D = np.zeros((3, nao, nao))
    Q = np.zeros((6, nao, nao))
    off = np.cumsum([0] + [2 * sh["l"] + 1 for sh in shells])
    for a, sha in enumerate(shells):
        for b, shb in enumerate(shells[:a + 1]):
            la, lb = sha["l"], shb["l"]
            nca = (la + 1) * (la + 2) // 2
            ncb = (lb + 1) * (lb + 2) // 2
            bS = np.zeros((nca, ncb))
            bD = np.zeros((3, nca, ncb))
            bQ = np.zeros((6, nca, ncb))
            for ca, aa in zip(sha["cn"], sha["exp"]):
                for cb, ab in zip(shb["cn"], shb["exp"]):
                    s0, d0, q0 = cart_moments(aa, ab, xyz[sha["at"]], xyz[shb["at"]],
                                              la, lb)
                    bS += ca * cb * s0
                    bD += ca * cb * d0
                    bQ += ca * cb * q0
            Ta, Tb = mats[la], mats[lb]
            sl_a = slice(off[a], off[a + 1])
            sl_b = slice(off[b], off[b + 1])
            sph = Ta @ bS @ Tb.T
            S[sl_a, sl_b] = sph
            S[sl_b, sl_a] = sph.T
            for d in range(3):
                m = Ta @ bD[d] @ Tb.T
                D[d, sl_a, sl_b] = m
                D[d, sl_b, sl_a] = m.T
            for qi in range(6):
                m = Ta @ bQ[qi] @ Tb.T
                Q[qi, sl_a, sl_b] = m
                Q[qi, sl_b, sl_a] = m.T
    return S, D, Q


def _pyscf_ref(zs, xyz_bohr, shells):
    from pyscf import gto
    SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
    atom, basis = [], {}
    for i, z in enumerate(zs):
        name = f"{SYM[z]}{i}"
        atom.append([name, tuple(xyz_bohr[i])])
        bs = [[sh["l"]] + [[float(e), float(c)] for e, c in zip(sh["exp"], sh["coef"])]
              for sh in shells if sh["at"] == i]
        basis[name] = bs
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", charge=0, spin=None, cart=False)
    r = mol.intor("int1e_r")                 # (3, nao, nao)
    rr = mol.intor("int1e_rr").reshape(3, 3, mol.nao, mol.nao)
    q6 = np.array([rr[a, b] for a, b in _Q6])
    return r, q6


def _gate():
    import math as _m
    BOHR = _K.BOHR
    ang = 104.5 * _m.pi / 180
    r_oh = 0.9572 * BOHR
    W = [[0.0, 0.0, 0.0],
         [r_oh * _m.sin(ang / 2), 0.0, r_oh * _m.cos(ang / 2)],
         [-r_oh * _m.sin(ang / 2), 0.0, r_oh * _m.cos(ang / 2)]]
    a4 = 1.087 * BOHR / _m.sqrt(3)
    cases = {
        "h2o": ([8, 1, 1], np.array(W)),
        "hf": ([1, 9], np.array([[0, 0, 0], [0, 0, 1.733]], float)),
        "ch4": ([6, 1, 1, 1, 1],
                np.array([[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4],
                          [-a4, a4, -a4], [-a4, -a4, a4]], float)),
    }
    fails = 0
    for name, (zs, xyz) in cases.items():
        shells, _ = OV.build_shells(zs, xyz, charge=0)
        S, D, Q = moment_matrices(zs, xyz, shells=shells, ao_order="pyscf")
        rref, qref = _pyscf_ref(zs, xyz, shells)
        dD = float(np.abs(D - rref).max())
        dQ = float(np.abs(Q - qref).max())
        ok = dD <= 1e-10 and dQ <= 1e-10
        fails += 0 if ok else 1
        print(f"  {name:5s} nao {S.shape[0]:2d}  max|dD| {dD:.2e}  max|dQ| {dQ:.2e}  "
              f"{'MATCH' if ok else 'MISMATCH'}")
    if fails:
        print(f"\nMULTIPOLE-INTEGRAL GATE: {fails} MISMATCH")
        sys.exit(1)
    print("\nMULTIPOLE-INTEGRAL GATE PASSED: D and Q match PySCF on the adapted basis")


if __name__ == "__main__":
    _gate()
