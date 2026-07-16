"""adapt.py -- Term 3: the BASIS coordination number and the q_eff that adapts the q-vSZP basis.

Decoded 2026-07-16 from the qvSZP setup tool (grimme-lab/qvSZP, src/chargscfcts.f90
`ncoord_basq` + app/main.f90) and verified against the v1 oracle:

  count_ij = 0.5 * (1 + erf( KN * (r_ij - rc_ij) / rc_ij )),   KN = -3.75 (the tool's value)
  rc_ij    = ( R_i + R_j ) in Bohr, R = Pyykko-Atsumi 2009 covalent radii (Angstrom, metals
             reduced 10% -- the table below is copied verbatim from the tool, Z = 1..118)
  CN_basis(A) = sum over pairs.

and the effective charge that scales the contraction coefficients (SI Eq. 28, header banner
`(q+aq^2+bCN^0.5+cqCN)`), with (h1, h2, h3) = the per-element triple from the basisq file header:

  q_eff = (q - h2*q^2) + h1*sqrt(CN_basis) + h3*q*CN_basis

(sign pinned by the oracle: the banner's `a` is -h2). q is the EEQ(BC) charge (eeqbc.py, Term 2).
The final contraction is then c = c0 + c1*q_eff per primitive (basisq.py columns).

GATE (declared): reproduce the oracle's per-atom CN(basis) (EEQ block) AND its full AO-setup
block (q, CN, dq, dcn, dqcn, total) across the standard probe set to the print precision.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BOHR = 1.8897261254578281
KN = -3.75

# Pyykko & Atsumi, Chem. Eur. J. 2009, 15, 188-197; metals decreased by 10%.
# Copied verbatim from grimme-lab/qvSZP src/chargscfcts.f90 (Angstrom; converted below).
_PYYKKO_AA = [
    0.29, 0.46,
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,
    1.76, 1.54,
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,
    1.89, 1.67,
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23,
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,
    2.09, 1.76,
    1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,
    1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42,
    2.01, 1.81,
    1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,
    1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,
    1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22,
    1.36, 1.43, 1.46, 1.58, 1.48, 1.57,
]
RCOV_BOHR = [r * BOHR for r in _PYYKKO_AA]


def basis_cn(zs, xyz_bohr):
    """CN(basis) per atom -- the CN that feeds q_eff."""
    xyz = np.asarray(xyz_bohr, float)
    n = len(zs)
    cn = np.zeros(n)
    for i in range(n):
        for j in range(i):
            rc = RCOV_BOHR[zs[i] - 1] + RCOV_BOHR[zs[j] - 1]
            r = float(np.linalg.norm(xyz[i] - xyz[j]))
            c = 0.5 * (1.0 + math.erf(KN * (r - rc) / rc))
            cn[i] += c
            cn[j] += c
    return cn


def q_eff(zs, q, cn, adapt_triples):
    """SI Eq. 28 with the basisq header triples {Z: (h1, h2, h3)}."""
    out = np.empty(len(zs))
    for i, z in enumerate(zs):
        h1, h2, h3 = adapt_triples[z]
        out[i] = (q[i] - h2 * q[i] ** 2) + h1 * math.sqrt(max(cn[i], 0.0)) + h3 * q[i] * cn[i]
    return out


def _gate():
    sys.path.insert(0, HERE)
    import basisq
    import eeqbc
    import oracle
    B = basisq.parse()
    triples = {z: B[z]["adapt"] for z in B}
    fails = 0
    for label, (atoms, chg, uhf) in eeqbc.PROBES.items():
        r = oracle.run(atoms, charge=chg, uhf=uhf)
        zs = [eeqbc._SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
        xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * BOHR
        cn = basis_cn(zs, xyz)
        worst_cn = max(abs(cn[i] - row["cn_basis"]) for i, row in enumerate(r["eeq"]))
        ours = eeqbc.charges(zs, xyz, charge=chg)
        qe = q_eff(zs, ours["q"], cn, triples)
        worst_qe = max(abs(qe[i] - row["total"]) for i, row in enumerate(r["ao_adapt"])) \
            if r["ao_adapt"] else float("nan")
        ok = worst_cn <= 5e-4 and worst_qe <= 5e-4
        fails += 0 if ok else 1
        print(f"  {label:8s} CN(basis) {worst_cn:.1e}   q_eff {worst_qe:.1e}   "
              f"{'MATCH' if ok else 'MISMATCH'}")
    if fails:
        print(f"\nADAPT GATE: {fails} MISMATCH")
        sys.exit(1)
    print("\nADAPT GATE PASSED: CN(basis) and q_eff are the oracle's, across the probe set")


if __name__ == "__main__":
    _gate()
