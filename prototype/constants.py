"""constants.py -- the ONE loader for every method constant that is not in the official files.

House rule (user-set, 2026-07-16): no method constant lives in code. Everything decoded from
reference sources or measured from the oracle sits in data/derived-constants.json with
provenance, is loaded ONCE here at import, and every prototype module reads it from this
namespace. When our Fortran port lands in src/tblite/, the same table becomes a parameter-file
extension -- our fork loads at startup what upstream hardcodes.

    import constants as K
    K.BASIS_CN_KN, K.BASIS_RCOV_BOHR[z], K.EEQBC["kbc"], K.REFOCC[z], K.D_PERM, ...
"""
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_PATH = os.path.join(_HERE, "data", "derived-constants.json")
_D = json.load(open(_PATH))

BOHR = 1.8897261254578281
EV = 27.211386245988          # Hartree -> eV (CODATA); the oracle prints eigenvalues in eV

# --- basis CN (q-vSZP) ----------------------------------------------------------------------
BASIS_CN_KN = _D["basis_cn"]["kn"]
BASIS_RCOV_BOHR = [r * BOHR for r in
                   _D["basis_cn"]["rcov_angstrom_pyykko_atsumi_2009_metals_minus_10pct"]]

# --- EEQ(BC) factory constants --------------------------------------------------------------
EEQBC = _D["eeqbc"]

# --- repulsion structural constants ---------------------------------------------------------
REPULSION_KEXP = _D["repulsion"]["kexp"]

# --- fractional reference occupations -------------------------------------------------------
_LMAP = {"s": 0, "p": 1, "d": 2, "f": 3}
REFOCC = {int(z): {_LMAP[l]: v for l, v in shells.items()}
          for z, shells in _D["reference_occupations"]["by_z"].items()}

# --- oracle conventions ---------------------------------------------------------------------
D_PERM = tuple(_D["oracle_conventions"]["d_perm_from_pyscf_m_order"])
D_SIGNS = tuple(_D["oracle_conventions"]["d_signs"])

# --- measured-but-not-closed onsite functions (calibration targets) -------------------------
ONSITE_MEASURED = _D["onsite_electronic_measured"]


def raw():
    """The full parsed table, for tooling."""
    return _D


if __name__ == "__main__":
    assert BASIS_CN_KN == -3.75 and len(BASIS_RCOV_BOHR) == 118
    assert abs(BASIS_RCOV_BOHR[0] - 0.29 * BOHR) < 1e-12
    assert EEQBC["kbc"] == 0.60 and EEQBC["norm_exp"] == 0.75
    assert REFOCC[8][0] == 1.673 and REFOCC[8][1] == 4.326
    assert REFOCC[14][2] == 0.287
    assert D_PERM == (2, 4, 1, 3, 0)
    print("CONSTANTS LOADER SELFTEST PASSED "
          f"({len(REFOCC)} refocc elements, {len(BASIS_RCOV_BOHR)} radii)")
