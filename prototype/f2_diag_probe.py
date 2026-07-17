"""f2_diag_probe.py -- the pz/px diagonal split: per-AO vs per-shell population resolution.

Measured at r_e: diag(pz) = -0.519, diag(px) = -0.579 -- split +0.060. Under a PER-SHELL
convention, level/mu/X/ES2 are identical for pz and px: the split must equal the ACP
diagonal difference alone. Under PER-AO, the exchange adds -gamma_pp*(m_pz - m_px) and the
ES2 shifts add their per-AO parts. Compute both predictions from the restart density and
the decoded kernels; the measured split decides.

Kernel guess for the cross-shell onsite exchange (to be refined): gamma_on(l,l') =
c_x(F)*(L5_l + L5_l') -- reduces to 2*c_x*L5 for l = l' (the atom-gated same-shell form).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import f2_stretch  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402

P_ = params.parse()
eF = P_["element"][9]
L5s, L5p = eF["shells"][3][0], eF["shells"][3][1]
CX = 0.6693 / 9.59                     # s(F)/9.59
G_SS = 2 * CX * L5s
G_SP = CX * (L5s + L5p)
G_PP = 2 * CX * L5p


def main():
    for R in (2.668, 3.8):
        rec = fock_recon.fock_ao([("F", 0, 0, 0), ("F", 0, 0, R / K.BOHR)])
        F, S = rec["F"], rec["S"]
        Ps = rec["state"]["P"] / 2.0
        m = np.diag(Ps @ S)                       # per-AO per-spin Mulliken populations
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        A = f2_stretch.acp_matrix([9, 9], xyz)
        meas_split = F[3, 3] - F[1, 1]
        acp_split = A[3, 3] - A[1, 1]
        # per-AO exchange difference: v_pz - v_px = gamma_pp*(m_pz - m_px)
        #   (+ gamma_sp*m_s identical for both -> cancels in the split)
        x_split_perAO = -G_PP * (m[3] - m[1])
        print(f"R = {R}:")
        print(f"  per-AO m: s {m[0]:.4f}  px {m[1]:.4f}  py {m[2]:.4f}  pz {m[3]:.4f}")
        print(f"  measured split (pz - px): {meas_split:+.5f}")
        print(f"  ACP split:                {acp_split:+.5f}")
        print(f"  per-shell prediction:     {acp_split:+.5f}  "
              f"(miss {meas_split - acp_split:+.5f})")
        print(f"  per-AO prediction:        {acp_split + x_split_perAO:+.5f}  "
              f"(X-part {x_split_perAO:+.5f}; miss "
              f"{meas_split - acp_split - x_split_perAO:+.5f})")
        # full diagonals, SAME-SHELL-ONLY exchange (the atom-gate form: Vx = -2c_x*L5_l*m)
        # + per-AO populations (the probe's verdict):
        mu_s, mu_p = -eF["shells"][5][0], -eF["shells"][5][1]
        lvl_s, lvl_p = -eF["shells"][0][0], -eF["shells"][0][1]
        for lab, i, lvl, mu, g in (("s ", 0, lvl_s, mu_s, G_SS),
                                   ("px", 1, lvl_p, mu_p, G_PP),
                                   ("pz", 3, lvl_p, mu_p, G_PP)):
            pred = lvl + A[i, i] + mu - g * m[i]
            print(f"  {lab}-diag: measured {F[i, i]:+.5f}  pred(lvl+ACP+mu+X_sameshell) "
                  f"{pred:+.5f}  miss {F[i, i] - pred:+.5f}")


if __name__ == "__main__":
    main()
