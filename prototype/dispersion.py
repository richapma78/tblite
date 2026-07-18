"""dispersion.py -- the g-xTB dispersion (gap #6), via the STOCK dftd4 library.

g-xTB uses self-consistent DFT-revD4 (SI Sec 1.17): a MODIFIED D4 with a sigmoidal charge
scaling zeta driven by the converged MULLIKEN charges (Eq 166-169), BJ damping with a2
eliminated (Eq 170), and Chai-Head-Gordon ATM. The four global damping params, DERIVED and
falsifiable from the parameter file (see derived-constants "dispersion_revd4_slots_DERIVED"):

    a1 = g1[9] = 1.2154627292      (BJ damping scale; a2 is eliminated -> a2 = 0)
    s8 = g2[9] = 0.304294728       (s6 = 1, s9 = 1 presumed; s9 ~ irrelevant at these sizes)

The REVERSED assignment misses by 100-2540 mEh, so this ordering is well-pinned.

WHAT THIS IS: stock dftd4 (charge+CN-dependent D4) with the g-xTB damping. It reproduces the
binary's printed revD4 to ~0.15-0.30 mEh on the six H..F references -- UNDER the 1 mEh
substitution grade, but NOT bit-exact: the residual is revD4's own sigmoidal Mulliken-charge
zeta (SI Eq 166-169), which stock D4 does not have. So this is the dispersion LOWER-BOUND /
stock approximation; closing the last ~0.3 mEh needs the zeta (a separate port step, needs the
atomic reference polarizabilities alpha+/alpha0/alpha-). The Fortran port will call dftd4
for real, so this is a faithful staging of that call.

Needs the dftd4 python bindings (pip install dftd4; present in gpudft and ~/dftd4env).
"""
import numpy as np

import paramfile

try:
    from dftd4.interface import DampingParam, DispersionModel
    _HAVE_D4 = True
except ImportError:                                    # keep the ledger runnable without D4
    _HAVE_D4 = False

_P = paramfile.load()
A1 = _P["g1"][9]                                       # 1.2154627292
S8 = _P["g2"][9]                                       # 0.304294728
S6, S9, A2 = 1.0, 1.0, 0.0                             # a2 eliminated in revD4 (SI Eq 170)


def available():
    return _HAVE_D4


def energy(zs, xyz, charge=0):
    """Stock-D4 dispersion (Eh) with the g-xTB damping. xyz in BOHR.
    Returns None if the dftd4 bindings are absent."""
    if not _HAVE_D4:
        return None
    model = DispersionModel(np.asarray(zs), np.asarray(xyz, float), charge=charge)
    res = model.get_dispersion(
        DampingParam(s6=S6, s8=S8, a1=A1, a2=A2, s9=S9), grad=False)
    return float(res["energy"])
