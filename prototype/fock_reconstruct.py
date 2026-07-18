"""fock_reconstruct.py -- the binary's CONVERGED total Fock, reconstructed EXACTLY from the
restart file, with no eigensolver extraction.

    F_ao = S . C^T . diag(eps) . C . S

where C = the restart MO coefficients (rows = MOs, orthonormal wrt S: C S C^T = I), eps = the
printed orbital energies (Eh), S = our (machine-exact) overlap. This holds because
F C_col = S C_col diag(eps) with C_col = C^T the column-MO matrix.

USE: F_ao is the binary's H0 + ACP + (all electron Fock terms) at the converged density. It is
the ground truth for the total Fock. To isolate the bare EHT Hamiltonian H0 you must still
subtract ACP (we have it, ~exact) AND the electron Fock in the binary's convention -- the
latter is the open piece: GE.fock's exchange potential is the old H2-grade skeleton, NOT the
validated mfx, so it does not subtract cleanly. Building the mfx exchange Fock (+ the ES1/ES2/
ES3 potentials as matrices) is what turns this into a clean H0 readout.

Even without that, F_ao is directly useful: its OFF-DIAGONALS are dominated by H0 (the electron
offdiagonals -- exchange, ES2 offsite -- are smaller), so F_ao offdiag vs GE H0 offdiag
localizes bonding-amplitude errors to ~10%. See PORT_STATUS.md / derived-constants
h0_electronic_gap_ROADMAP.
"""
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import gxtb_engine as GE  # noqa: E402
import restart  # noqa: E402

BOHR = 1.8897261254578281


def binary_fock(zs, xyz_bohr, sym):
    """F_ao (Eh) at the converged density, plus S/H0/A/meta from GE. xyz in bohr."""
    atoms = [(sym[z], x / BOHR, y / BOHR, z2 / BOHR) for z, (x, y, z2) in zip(zs, xyz_bohr)]
    st = restart.converged_state(atoms)
    n = st["nsao"]
    C = st["C"].reshape(n, n)
    eps = np.array(sorted(st["eps_ev"])) / 27.211386245988      # eV -> Eh, ascending
    B = GE.build(zs, np.array(xyz_bohr, float))
    S = B["S"]
    Fao = S @ C.T @ np.diag(eps) @ C @ S
    return {"Fao": Fao, "S": S, "H0": B["H0"], "A": B["A"], "meta": B["meta"],
            "eps": eps, "P": st["P"]}


if __name__ == "__main__":
    import os
    os.environ.setdefault("GXTB_H0V2", "1")
    np.set_printoptions(precision=4, suppress=True, linewidth=140)
    r = binary_fock([1, 9], np.array([[0, 0, 0], [0, 0, 1.733]], float), {1: "H", 9: "F"})
    print("AO order 0=H_s 1=F_s 2=F_px 3=F_py 4=F_pz")
    print("F_ao (binary converged Fock):\n", r["Fao"])
    print("GE H0_v2 + ACP:\n", r["H0"] + r["A"])
    print("H_s-F_s   F_ao %+.4f  GE %+.4f" % (r["Fao"][0, 1], (r["H0"] + r["A"])[0, 1]))
