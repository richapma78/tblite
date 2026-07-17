"""camm.py -- cumulative atomic multipole moments (AES step 2), validated against the
binary's printed dipoles. Two checkpoints g-xTB hands us for free:
  full molecular dipole  = Sum_A Zval_A R_A  -  Tr(P D)      (exact from the density)
  point-charge dipole    = Sum_A q_A R_A      (Mulliken charges only)
For H2O the binary prints 0.9289 au and 0.3107 au; matching both validates the D
integrals, the Mulliken charges, and isolates the atomic-dipole (CAMM) contribution --
the AES ingredient. Then the per-atom CAMM q_A / mu_A / theta_A are assembled (Mulliken
partition, shifted to atomic centres) for the kernel step."""
import math
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import multipole as MP  # noqa: E402
import overlap as OV  # noqa: E402
import restart  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
# integer valence-electron count = the effective nuclear charge for the electrostatics
# (NOT the fractional reference occupation -- that shift is the F/HF 0.24% dipole error)
ZVAL = {1: 1, 6: 4, 7: 5, 8: 6, 9: 7}


def systems():
    ang = 104.5 * math.pi / 180
    r_oh = 0.9572 * BOHR
    W = [[0.0, 0.0, 0.0],
         [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
         [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
    a4 = 1.087 * BOHR / math.sqrt(3)
    r_nh = 1.012 * BOHR
    st, ct = 0.9262, -0.3770
    return [
        ("H2O", [8, 1, 1], np.array(W)),
        ("HF", [1, 9], np.array([[0, 0, 0], [0, 0, 1.733]], float)),
        ("CH4", [6, 1, 1, 1, 1], np.array(
            [[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4], [-a4, a4, -a4], [-a4, -a4, a4]],
            float)),
        ("NH3", [7, 1, 1, 1], np.array(
            [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                            r_nh * st * math.sin(2 * math.pi * k / 3),
                            r_nh * ct] for k in range(3)], float)),
    ]


def printed_dipoles(raw):
    """the au dipole magnitudes g-xTB prints on the 'total (au/Debye)' lines."""
    return [float(m) for m in
            re.findall(r"total\s+\(au/Debye\)\s*:\s*(-?\d+\.\d+)", raw)]


for name, zs, xyz in systems():
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    st = restart.converged_state(atoms)
    P = st["P"]
    raw = st["raw"]
    shells, _ = OV.build_shells(zs, xyz, charge=0)
    S, D, Q = MP.moment_matrices(zs, xyz, shells=shells, ao_order="oracle")
    n = S.shape[0]
    # AO -> atom map (shell order matches build order)
    ao_at = []
    for sh in shells:
        ao_at += [sh["at"]] * (2 * sh["l"] + 1)
    ao_at = np.array(ao_at)
    # full molecular dipole from the density
    mu_elec = -np.einsum("xij,ij->x", D, P)
    mu_nuc = sum(ZVAL[zs[a]] * xyz[a] for a in range(len(zs)))
    mu_full = mu_nuc + mu_elec
    # Mulliken atomic charges -> point-charge dipole
    PS = P @ S
    pop = np.array([PS[i, i] for i in range(n)])
    qat = np.array([ZVAL[zs[a]] - pop[ao_at == a].sum() for a in range(len(zs))])
    mu_pt = sum(qat[a] * xyz[a] for a in range(len(zs)))
    got = printed_dipoles(raw)
    print(f"\n{name}: printed dipoles (au) = {got}")
    print(f"  full  |mu| ours {np.linalg.norm(mu_full):.5f}  vec {mu_full.round(4)}")
    print(f"  point |mu| ours {np.linalg.norm(mu_pt):.5f}  vec {mu_pt.round(4)}")
    print(f"  sum q = {qat.sum():+.6f}  (neutral check)   q_at {qat.round(4)}")
