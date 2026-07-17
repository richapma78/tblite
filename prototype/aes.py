"""aes.py -- the AES (anisotropic electrostatics) energy, gap #5, assembled from:
  - CAMM atomic moments q/dp/qp from the g-xTB density (validated via the dipole),
  - the tblite tensor structure (get_multipole_matrix_0d + get_energy_aes),
  - g-xTB's OWN erf-damped kernels with parameters extracted live from the binary.

Kernel (from mmomgaberf, gdb-extracted): gamma_n(R) = s_n * 0.5(1+erf(a_n(R-R0)))/R^n,
  s3=G2[5]=0.4667, s5=G2[6]=0.2287, s7=G2[7]=0.0523, s9=-0.003; a3=0.45, a5=a7=a9=1.05;
  R0 = per-element-pair table (bohr). tblite tensor forms: amat_sd = vec*K3;
  amat_dd = unity*K3d - vec@vec*3*K5; amat_sq = [xx,2xy,yy,2xz,2yz,zz]*K5.
  Energy: 0.5*e01(charge-dip) + e11(dip-dip) + 0.5*e02(charge-quad).

4-POINT GATE (declared): reproduce the binary's printed ES multipole on h2o/hf/ch4/nh3
(0.00405221 / 0.000011 / 0.002332 / 0.006833) -- a wrong assembly cannot match all four.
Convention flags (DD_DAMP, QUAD) are swept if the first mapping misses; the gate decides.
"""
import math
import re
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import multipole as MP  # noqa: E402
import overlap as OV  # noqa: E402
import restart  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
ZVAL = {1: 1, 6: 4, 7: 5, 8: 6, 9: 7}
# gdb-extracted AES parameters (pass 89)
A3, A5 = 0.45, 1.05
S3, S5, S7, S9 = 0.4667106684, 0.2287119112, 0.0522538831, -0.003
R0 = {(1, 1): 2.1823, (1, 6): 2.4492, (1, 7): 2.3667, (1, 8): 2.1768, (1, 9): 2.0646,
      (6, 6): 2.9103, (6, 7): 2.7063, (6, 8): 2.5697, (6, 9): 2.4770, (7, 7): 2.6225,
      (7, 8): 2.4846, (7, 9): 2.3885, (8, 8): 2.4817, (8, 9): 2.3511, (9, 9): 2.2996}
PRINTED = {"h2o": 0.00405221, "hf": 0.000011, "ch4": 0.002332, "nh3": 0.006833}


def r0_of(zi, zj):
    return R0[(min(zi, zj), max(zi, zj))]


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
        "h2o": ([8, 1, 1], np.array(W)),
        "hf": ([1, 9], np.array([[0, 0, 0], [0, 0, 1.733]], float)),
        "ch4": ([6, 1, 1, 1, 1], np.array(
            [[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4], [-a4, a4, -a4], [-a4, -a4, a4]],
            float)),
        "nh3": ([7, 1, 1, 1], np.array(
            [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                            r_nh * st * math.sin(2 * math.pi * k / 3),
                            r_nh * ct] for k in range(3)], float)),
    }


# quadrupole packing order [xx, xy, yy, xz, yz, zz] (tblite), off-diag doubled in amat_sq
_QORD = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]
_MPSCALE = np.array([1., 2., 1., 2., 2., 1.])  # off-diagonals count twice in contraction


def camm(zs, xyz, P, S, D, Q, ao_at):
    """per-atom CAMM (SI Eq 112/113b): qat(nat), dpat(3,nat), Theta(nat,3,3) traceless."""
    nat = len(zs)
    n = S.shape[0]
    PS = P @ S
    pop = np.array([PS[i, i] for i in range(n)])
    qat = np.array([ZVAL[zs[a]] - pop[ao_at == a].sum() for a in range(nat)])
    dpat = np.zeros((3, nat))
    Theta = np.zeros((nat, 3, 3))
    q_src = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (1, 1): 3, (1, 2): 4, (2, 2): 5}
    for a in range(nat):
        mu_on = np.where(ao_at == a)[0]
        Ra = xyz[a]
        Pa = P[mu_on, :]
        for al in range(3):                                     # Eq 112b
            dpat[al, a] = np.sum(Pa * (Ra[al] * S[mu_on, :] - D[al, mu_on, :]))
        raw = np.zeros((3, 3))
        for (a1, a2) in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]:
            qs = q_src[(a1, a2)]                                 # Eq 112c
            val = np.sum(Pa * (Ra[a1] * D[a2, mu_on, :] + Ra[a2] * D[a1, mu_on, :]
                               - Ra[a1] * Ra[a2] * S[mu_on, :] - Q[qs, mu_on, :]))
            raw[a1, a2] = raw[a2, a1] = val
        Theta[a] = 1.5 * raw - 0.5 * np.trace(raw) * np.eye(3)   # Eq 113b traceless
    # SI Eq-116 sign convention for the multipole moments (physical, electron charge
    # negative): the CAMM density moments carry the opposite sign of Eq 112b's Mulliken
    # form. Flips the odd-order (charge-dipole, dipole-quad) terms; the gate confirms.
    return qat, -dpat, Theta


def energy(zs, xyz, qat, dpat, Theta):
    """SI Eq 116: the full damped multipole energy up to quad-quad, with the
    gdb-extracted erf kernels gamma_n = s_n*0.5(1+erf(a_n(R-R0)))/R^n."""
    nat = len(zs)
    E = 0.0
    for A in range(nat):
        for B in range(nat):
            if A == B:
                continue
            R = xyz[A] - xyz[B]
            r = float(np.linalg.norm(R))
            r2 = r * r
            r0 = r0_of(zs[A], zs[B])
            g3 = S3 * 0.5 * (1 + erf(A3 * (r - r0))) / r ** 3
            g5 = S5 * 0.5 * (1 + erf(A5 * (r - r0))) / r ** 5
            g7 = S7 * 0.5 * (1 + erf(A5 * (r - r0))) / r ** 7
            g9 = S9 * 0.5 * (1 + erf(A5 * (r - r0))) / r ** 9
            qA, qB = qat[A], qat[B]
            mA, mB = dpat[:, A], dpat[:, B]
            TA, TB = Theta[A], Theta[B]
            mAR, mBR = mA @ R, mB @ R
            RTA = R @ TA @ R
            RTB = R @ TB @ R
            # charge-dipole (f3)
            E += 0.5 * g3 * (mAR * qB - qA * mBR)
            # dipole-dipole (f5)
            E += -0.5 * g5 * (3 * mAR * mBR - (mA @ mB) * r2)
            # charge-quadrupole (f5)
            E += 0.5 * g5 * (qB * RTA + qA * RTB)
            # dipole-quadrupole (f7)
            E += -0.5 * g7 * (5 * RTA * mBR - 2 * (R @ TA @ mB) * r2
                              + 5 * mAR * RTB - 2 * (R @ TB @ mA) * r2)
            # quadrupole-quadrupole (f9)
            E += (1.0 / 6.0) * g9 * (35 * RTA * RTB - 20 * (R @ TA @ TB @ R) * r2
                                     + 3 * np.sum(TA * TB) * r2 * r2)
    return float(E)


def run(verbose=True):
    worst = 0.0
    for name, (zs, xyz) in systems().items():
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        st = restart.converged_state(atoms)
        P = st["P"]
        # printed ES multipole from the SAME run (exact, geometry-consistent)
        printed = float(re.search(r"ES multipole\s*:\s*(-?\d+\.\d+)",
                                  st["raw"]).group(1))
        shells, _ = OV.build_shells(zs, xyz, charge=0)
        S, D, Q = MP.moment_matrices(zs, xyz, shells=shells, ao_order="oracle")
        ao_at = np.array([sh["at"] for sh in shells
                          for _ in range(2 * sh["l"] + 1)])
        qat, dpat, Theta = camm(zs, xyz, P, S, D, Q, ao_at)
        e = energy(zs, xyz, qat, dpat, Theta)
        d = e - printed
        PRINTED[name] = printed
        worst = max(worst, abs(d))
        if verbose:
            print(f"  {name:4s} ours {e:+.6f}  printed {PRINTED[name]:+.6f}  "
                  f"d {d:+.6f}")
    if verbose:
        print(f"  worst |d| = {worst:.6f}  "
              f"({'PASS' if worst <= 5e-5 else 'MISS'} bar 5e-5)")
    return worst


if __name__ == "__main__":
    run()
