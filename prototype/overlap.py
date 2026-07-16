"""overlap.py -- Term 4: overlap integrals over the ADAPTED q-vSZP basis, gated against PySCF.

The basis for a given molecule is built by the gated chain: EEQ(BC) charges (eeqbc.py) ->
CN(basis) + q_eff (adapt.py) -> per-primitive contraction c = c0 + c1*q_eff (basisq.py columns).
This module evaluates the overlap matrix over that basis with our own integral engine:

  - Cartesian primitive overlaps by the Obara-Saika 1-D recursion (exact, all angular momenta)
  - real-spherical transformation built PROGRAMMATICALLY from the solid-harmonic closed form
    (no hand-copied tables to mistype), in PySCF's m-ordering
  - PySCF's normalization convention mirrored exactly (coefficients refer to normalized
    primitives; each contracted shell normalized to unit self-overlap), so the external gate
    compares like with like.

GATE (declared): our S matrix equals PySCF's `int1e_ovlp` for the SAME adapted basis fed to a
per-atom custom-basis Mole, max |dS| <= 1e-10, on probes covering s/p (water), heavier p + d
(PdCl2), and an f-shell species if the basis carries one. The ORACLE's normalization convention
is deliberately NOT asserted here -- it gets pinned by the eigenvalue/population gates once H0
exists; this gate proves the INTEGRALS.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# ---------------------------------------------------------------- cartesian primitive overlap
def _os_1d(la, lb, pa, pb, mu):
    """Obara-Saika 1-D overlap table s[i,j] (unnormalized, without the Gaussian prefactor)."""
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


def cart_prims(a, b, ra, rb, la, lb):
    """Overlap block over CARTESIAN components for one primitive pair (unnormalized prims)."""
    mu = a + b
    p = (a * ra + b * rb) / mu
    ab2 = float(np.dot(ra - rb, ra - rb))
    pref = math.exp(-a * b / mu * ab2) * (math.pi / mu) ** 1.5
    tabs = [_os_1d(la, lb, p[k] - ra[k], p[k] - rb[k], mu) for k in range(3)]
    ca, cb = cart_components(la), cart_components(lb)
    out = np.empty((len(ca), len(cb)))
    for i, (ax, ay, az) in enumerate(ca):
        for j, (bx, by, bz) in enumerate(cb):
            out[i, j] = pref * tabs[0][ax, bx] * tabs[1][ay, by] * tabs[2][az, bz]
    return out


def cart_components(l):
    """PySCF's cartesian component ordering: lexicographic with x >= y >= z pattern."""
    return [(l - i, i - j, j) for i in range(l + 1) for j in range(i + 1)]


# ---------------------------------------------------------- real-spherical transform (general)
def _c2s(l):
    """Cartesian->real-spherical matrix (2l+1, ncart), rows m = -l..+l (PySCF's ordering).

    Built NUMERICALLY, not from a transcribed table: the real spherical harmonic Y_lm restricted
    to the unit sphere is an exact linear combination of the degree-l monomials (the monomial
    space on the sphere has dimension (l+1)(l+2)/2 = (2l+1) + (2l-3) + ... so the representation
    is unique), so a least-squares solve on random unit vectors recovers the coefficients to
    machine precision. Scaling: each row is normalized so that a shell of AXIS-NORMALIZED
    cartesian primitives (prim_norm below -- PySCF's gto_norm convention) yields unit-norm
    spherical functions; the (2l-1)!!-vs-component factors are absorbed here, which is exactly
    the CCA/PySCF convention the gate compares against."""
    try:
        from scipy.special import sph_harm_y

        def _ylm(m, l_, phi, theta):
            return sph_harm_y(l_, m, theta, phi)        # new scipy: (n, m, polar, azimuth)
    except ImportError:
        from scipy.special import sph_harm

        def _ylm(m, l_, phi, theta):
            return sph_harm(m, l_, phi, theta)          # old scipy: (m, n, azimuth, polar)
    rng = np.random.default_rng(l + 7)
    npts = 4 * (l + 1) * (l + 2)
    u = rng.normal(size=(npts, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    theta = np.arccos(np.clip(u[:, 2], -1, 1))          # polar
    phi = np.arctan2(u[:, 1], u[:, 0])                  # azimuth
    comps = cart_components(l)
    M = np.stack([u[:, 0] ** lx * u[:, 1] ** ly * u[:, 2] ** lz for lx, ly, lz in comps], axis=1)
    rows = []
    for m in range(-l, l + 1):
        ylm = _ylm(abs(m), l, phi, theta)
        if m > 0:
            y = math.sqrt(2.0) * (-1) ** m * ylm.real
        elif m < 0:
            y = math.sqrt(2.0) * (-1) ** m * ylm.imag
        else:
            y = ylm.real
        coef, *_ = np.linalg.lstsq(M, y, rcond=None)
        rows.append(coef)
    T = np.array(rows)
    if l == 1:
        T = T[[2, 0, 1]]        # PySCF quirk: p shells ordered (x, y, z), not m = -1..+1
    # row scaling: with axis-normalized cartesian overlaps Scc (single primitive, exponent 1),
    # require diag(T Scc T^T) = 1
    a = 1.0
    n = np.array([prim_norm(a, l)] * len(comps))
    Scc = cart_prims(a, a, np.zeros(3), np.zeros(3), l, l) * np.outer(n, n)
    d = np.sqrt(np.diag(T @ Scc @ T.T))
    return T / d[:, None]


# -------------------------------------------------------------------------- shells and matrix
def prim_norm(a, l):
    """Norm of a cartesian primitive with all angular momentum on one axis (the shell norm
    PySCF uses via gto_norm): ((2a/pi)^(3/4)) * (4a)^(l/2) / sqrt((2l-1)!!)."""
    dfac = 1.0
    for k in range(2 * l - 1, 0, -2):
        dfac *= k
    return (2.0 * a / math.pi) ** 0.75 * (4.0 * a) ** (l / 2.0) / math.sqrt(dfac)


def build_shells(zs, xyz_bohr, charge=0):
    """The full gated chain -> per-shell (l, center, exps, ADAPTED contracted coeffs)."""
    import adapt
    import basisq
    import eeqbc
    B = basisq.parse()
    q = eeqbc.charges(zs, xyz_bohr, charge=charge)["q"]
    cn = adapt.basis_cn(zs, xyz_bohr)
    qe = adapt.q_eff(zs, q, cn, {z: B[z]["adapt"] for z in B})
    shells = []
    for i, z in enumerate(zs):
        for l, prims in B[z]["shells"]:
            exps = np.array([p[0] for p in prims])
            coef = np.array([p[1] + p[2] * qe[i] for p in prims])
            shells.append({"l": l, "at": i, "exp": exps, "coef": coef})
    return shells, qe


def overlap(zs, xyz_bohr, charge=0, shells=None):
    """Spherical-AO overlap matrix in PySCF ordering/normalization for the adapted basis."""
    xyz = np.asarray(xyz_bohr, float)
    if shells is None:
        shells, _ = build_shells(zs, xyz, charge=charge)
    mats = {l: _c2s(l) for l in sorted({sh["l"] for sh in shells})}
    # normalized-primitive coefficients, then contracted self-overlap normalization
    for sh in shells:
        c = sh["coef"] * np.array([prim_norm(a, sh["l"]) for a in sh["exp"]])
        l, e = sh["l"], sh["exp"]
        dfac = 1.0
        for k in range(2 * l - 1, 0, -2):
            dfac *= k
        ee = e[:, None] + e[None, :]
        s_rad = (np.pi ** 1.5) * dfac / (2.0 ** l) / ee ** (l + 1.5)
        norm2 = float(c @ s_rad @ c)
        sh["cn"] = c / math.sqrt(norm2)
    nao = sum(2 * sh["l"] + 1 for sh in shells)
    S = np.zeros((nao, nao))
    off = np.cumsum([0] + [2 * sh["l"] + 1 for sh in shells])
    for a, sha in enumerate(shells):
        for b, shb in enumerate(shells[:a + 1]):
            la, lb = sha["l"], shb["l"]
            blk = np.zeros(((la + 1) * (la + 2) // 2, (lb + 1) * (lb + 2) // 2))
            for ca, aa in zip(sha["cn"], sha["exp"]):
                for cb, ab in zip(shb["cn"], shb["exp"]):
                    blk += ca * cb * cart_prims(aa, ab, xyz[sha["at"]], xyz[shb["at"]], la, lb)
            sph = mats[la] @ blk @ mats[lb].T
            S[off[a]:off[a + 1], off[b]:off[b + 1]] = sph
            S[off[b]:off[b + 1], off[a]:off[a + 1]] = sph.T
    return S


# ------------------------------------------------------------------------------- PySCF gate
def _pyscf_S(zs, xyz_bohr, shells):
    from pyscf import gto
    SYM = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S", 17: "Cl", 46: "Pd", 58: "Ce"}
    atom, basis = [], {}
    for i, z in enumerate(zs):
        name = f"{SYM[z]}{i}"
        atom.append([name, tuple(xyz_bohr[i])])
        bs = []
        for sh in shells:
            if sh["at"] == i:
                bs.append([sh["l"]] + [[float(e), float(c)]
                                       for e, c in zip(sh["exp"], sh["coef"])])
        basis[name] = bs
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", charge=0, spin=None, cart=False)
    return mol.intor("int1e_ovlp")


def _gate():
    import eeqbc
    BOHRC = 1.8897261254578281
    cases = {k: eeqbc.PROBES[k] for k in ("water", "AcCl", "PdCl2")}
    # f-shell coverage: cerium carries an f(7) shell in this basis (lmax 3)
    cases["CeO"] = ([("Ce", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.82)], 0, 0)
    eeqbc._SYM2Z.setdefault("ce", 58)
    fails = 0
    for label, (atoms, chg, uhf) in cases.items():
        zs = [eeqbc._SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
        xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * BOHRC
        shells, _ = build_shells(zs, xyz, charge=chg)
        S = overlap(zs, xyz, charge=chg, shells=shells)
        Sref = _pyscf_S(zs, xyz, shells)
        d = float(np.abs(S - Sref).max())
        lmax = max(sh["l"] for sh in shells)
        ok = d <= 1e-10
        fails += 0 if ok else 1
        print(f"  {label:8s} nao {S.shape[0]:3d}  lmax {lmax}  max|dS| {d:.2e}  "
              f"{'MATCH' if ok else 'MISMATCH'}")
    if fails:
        print(f"\nOVERLAP GATE: {fails} MISMATCH")
        sys.exit(1)
    print("\nOVERLAP GATE PASSED: our integral engine equals PySCF on the adapted basis")


if __name__ == "__main__":
    _gate()
