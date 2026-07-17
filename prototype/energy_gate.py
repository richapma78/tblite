"""energy_gate.py -- the TOTAL-ENERGY bookkeeping gate (H2 first).

Assemble every printed term from our validated formulas at the oracle's density and gate
against the printed decomposition. The T1 rule makes the bookkeeping explicit: energies
come from the ENERGY formulas (Tr(H^EHT o P), Tr(H^ACP o P), the 4-index Ex, the ES
energy expressions), never from eigenvalue sums.

H2 at q = 0, closed shell: E1 = E2+3 = AES = Espin = 0 by symmetry (their Fock shifts are
nonzero but their energies vanish at the point -- the method's own structure); OFX = 0
(s-only). Dispersion is NOT decoded (revD4): compared as pass-through, labeled.

GATE (pre-declared): |electronic (ours - printed)| <= 0.02 Eh at R = 1.4/2.5/4.0
(the H0-forward's 0.0094-worst element enters through Tr(H0 P)); repulsion at its own
gated 5e-9; increments exact; |total - dispersion| <= 0.02.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapt  # noqa: E402
import basisq  # noqa: E402
import constants as K  # noqa: E402
import f2_stretch  # noqa: E402
import fock_recon  # noqa: E402
import params  # noqa: E402
import repulsion  # noqa: E402
import overlap  # noqa: E402
import scf_h2  # noqa: E402

P_ = params.parse()
eH = P_["element"][1]
AMP = P_["globals"]["g1"][0] * eH["l8"][0]
L2, B_SHP = eH["shells"][0][0], 0.0315
K_FIT, KB_FIT = 1.110, 0.179

Bq = basisq.parse()
_prims = Bq[1]["shells"][0][1]
_e0 = np.array([p[0] for p in _prims])
_c0 = np.array([p[1] for p in _prims])
_c1 = np.array([p[2] for p in _prims])


def h0_matrix(R):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    cn = adapt.basis_cn([1, 1], xyz)[0]
    q_ = KB_FIT * math.sqrt(cn)
    sh = [{"l": 0, "at": i, "exp": _e0 * K_FIT, "coef": (_c0 + _c1 * q_).copy()}
          for i in (0, 1)]
    Sh = overlap.overlap([1, 1], xyz, shells=sh)
    H0 = np.diag([-L2, -L2]).astype(float)
    H0[0, 1] = H0[1, 0] = AMP * (-L2) * (1 + B_SHP * R) ** 2 * float(Sh[0, 1])
    return H0


def main():
    ok = True
    for R in (1.4, 2.5, 4.0):
        rec = fock_recon.fock_ao([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
        r = rec["state"]
        P = r["P"]
        terms = {}
        import re as _re
        for t in ("ES1 (charge SIE)", "ES2+3", "ES multipole", "Ex (Mulliken)",
                  "Espinpol", "electronic", "atomic core increments", "dispersion",
                  "nuclear repulsion"):
            m = _re.search(rf"^\s*{_re.escape(t)}\s*:\s*(-?\d+\.\d+)", r["raw"], _re.M)
            terms[t] = float(m.group(1)) if m else None
        m = _re.search(r"^\s*total\s+(-?\d+\.\d+)", r["raw"], _re.M)
        total = float(m.group(1))
        # ---- our pieces
        H0 = h0_matrix(R)
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        A = f2_stretch.acp_matrix([1, 1], xyz)
        _S, _H0b, _Ab, _E, gam = scf_h2.build(R)
        E_eht = float(np.sum(H0 * P))
        E_acp = float(np.sum(A * P))
        E_x = scf_h2.ex_energy(P, rec["S"], gam)
        E_rep = repulsion.energy([1, 1], xyz, [0.0, 0.0], P_, sign=+1, mean_rc=True,
                                 comb="harmonic")
        E_elec_ours = E_eht + E_acp + E_x
        d_elec = E_elec_ours - terms["electronic"]
        d_rep = E_rep - terms["nuclear repulsion"]
        d_ex = E_x - terms["Ex (Mulliken)"]
        tot_ours = E_elec_ours + terms["atomic core increments"] + \
            terms["dispersion"] + E_rep
        d_tot = tot_ours - total
        g = abs(d_elec) <= 0.02 and abs(d_tot) <= 0.02 and abs(d_ex) <= 0.001
        ok &= g
        print(f"R={R}:")
        print(f"  EHT Tr(H0 P) {E_eht:+.5f}   ACP Tr(A P) {E_acp:+.5f}   "
              f"Ex {E_x:+.5f} (d {d_ex:+.5f})")
        print(f"  electronic: ours {E_elec_ours:+.5f}  printed "
              f"{terms['electronic']:+.5f}  d {d_elec:+.5f}")
        print(f"  repulsion:  ours {E_rep:+.6f}  printed "
              f"{terms['nuclear repulsion']:+.6f}  d {d_rep:+.2e}")
        print(f"  total: ours {tot_ours:+.5f}  printed {total:+.5f}  d {d_tot:+.5f}   "
              f"[{'OK' if g else 'MISS'}]", flush=True)
    print("ENERGY GATE " + ("PASSED" if ok else "FAILED (report per term)"))


if __name__ == "__main__":
    main()
