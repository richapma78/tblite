"""es2_energy.py -- the DECODED second-order electrostatic energy (SI Eq 100-102),
replacing the interim es_charge_energy in the ledger.

  E2   = 0.5 sum_l sum_l' q_l q_l' gamma2_ll'
  gamma2_lAlB = 1 / [ R_AB + 0.5 (1/U_lA + 1/U_lB) exp(-k2x R_AB) ]      (Eq 101)
  U_lA        = T32_lA * ipse_A * (1 + Gamma_A * CN_A)                    (Eq 102)
  CN_A        = sqrt( sum_B p_AB^2 ),  p = 0.5(1+erf(-kcn (R-rc)/rc))     (internal L2/sqrt CN)

Every method constant is sourced from gxtb_parameters via paramfile; k2x=g2[1], kcn=g1[5].
rcov is the derived covalent-radius table (see internal_CN_is_L2norm_sqrtCN_DERIVED; flagged
for exact Ghidra sourcing). Verified: with the binary's OWN gamma2 the sum is bit-exact
(es2_gamma_gate); with the computed CN it reproduces the binary E2 to ~sub-0.1 mEh.
"""
import math

import paramfile

_P = paramfile.load()
KCN = _P["g1"][5]
K2X = _P["g2"][1]
# derived covalent radii (bohr) -- worst |dp|=2.3e-3 over the HF/FF sweeps + molecules
RCOV = {1: 0.6650, 6: 1.1465, 7: 1.0264, 8: 0.9288, 9: 0.7815}


def _nsh(z):
    return sum(1 for t in _P["elements"][z]["T32"] if t != 0.0)


def coordination(zs, xyz):
    """Internal CN_A = sqrt(sum_B p_AB^2), p the erf counting function."""
    n = len(zs)
    cn = [0.0] * n
    for i in range(n):
        s = 0.0
        for j in range(n):
            if i == j:
                continue
            R = math.dist(xyz[i], xyz[j])
            rc = RCOV[zs[i]] + RCOV[zs[j]]
            s += (0.5 * (1.0 + math.erf(-KCN * (R - rc) / rc))) ** 2
        cn[i] = math.sqrt(s)
    return cn


def _U(z, l, cn_a):
    e = _P["elements"][z]
    return e["T32"][l] * e["ipse"] * (1.0 + e["gamma"] * cn_a)


def energy(zs, xyz, q, cn=None):
    """E2 in Eh. zs: atomic numbers; xyz: (n,3) bohr; q: {(atom, l): shell charge}.
    cn: optional per-atom CN (else computed). Pass the binary's CN to isolate the sum."""
    if cn is None:
        cn = coordination(zs, xyz)
    shells = [(a, l, z) for a, z in enumerate(zs) for l in range(_nsh(z))]
    E = 0.0
    for (ai, li, zi) in shells:
        Ui = _U(zi, li, cn[ai])
        qi = q[(ai, li)]
        for (aj, lj, zj) in shells:
            Uj = _U(zj, lj, cn[aj])
            R = math.dist(xyz[ai], xyz[aj]) if ai != aj else 0.0
            g = 1.0 / (R + 0.5 * (1.0 / Ui + 1.0 / Uj) * math.exp(-K2X * R))
            E += 0.5 * qi * q[(aj, lj)] * g
    return E


if __name__ == "__main__":
    # HF at 1.733 bohr, binary rdx charges -> isolate gamma2/CN accuracy vs E2=0.04502243
    zs = [1, 9]
    xyz = [(0.0, 0.0, 0.0), (0.0, 0.0, 1.733)]
    q = {(0, 0): 0.412045, (1, 0): -0.093354, (1, 1): -0.318692}
    cn = coordination(zs, xyz)
    print(f"computed CN (HF) = {cn[0]:.5f}, {cn[1]:.5f}  (binary 0.30617)")
    e_comp = energy(zs, xyz, q)
    e_bincn = energy(zs, xyz, q, cn=[0.3061652499, 0.3061652499])
    print(f"E2 (computed CN)   = {e_comp:.8f}  d vs binary {(e_comp-0.04502243)*1000:+.4f} mEh")
    print(f"E2 (binary CN)     = {e_bincn:.8f}  d vs binary {(e_bincn-0.04502243)*1000:+.4f} mEh")
    print(f"binary E2 (rdx)    = 0.04502243   printed ES2+3 = 0.04504671")
