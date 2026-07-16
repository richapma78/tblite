"""params.py -- the gxtb_parameters file, parsed with NAMES where the perturbation map pinned them.

Layout (structural census + perturbation mapping, 2026-07-16; data/perturb_map.json is the
evidence): 2 x 10 GLOBALS, then per-element blocks "Z" + rows L1(10) L2..L7(4 each, per-shell
s/p/d/f) L8(8) L9(8). 79 elements parameterized (H..U with 4f/5f gaps: missing 59-69, 90-91,
93-103).

Named so far (slot -> role; evidence = which oracle observable moved when ONLY that slot did):
  L1[0] alpha0    repulsion exponent      (rep-only; far-range response grows like (R+R0)^1.5)
  L1[1] r0        repulsion offset radius (rep-only; far-range response grows like (R+R0)^0.5)
  L1[2] zeff0     effective nuclear charge (rep-only, positive response, grows with Z)
  L1[3] kcn_rep   CN-dependence of alpha  (rep-only, DEAD at CN~0, strong in NH4+)
  L1[4] kq_rep    charge-dependence of Zeff (rep-only, dead for q=0 pairs)
  L1[5] rcov_cn   Eq.47 covalent radius (the INTERNAL/Hamiltonian CN: moves electronic terms
                  AND repulsion-via-alpha(CN) together)
  L1[6] ku_cn    Hubbard CN-dependence U = U0*(1+ku*CN) (SI Eq. 102; shape-identified on the
                  HCl stretch: small ES23+ES1 response, erf-fast decay, dead by R=5)
  L1[7] hbasis   Hamiltonian-basis scale (persists at dissociation BY MOVING THE CHARGES --
                  the one L1 slot whose response never dies with R; Ex-heavy at short range)
  L1[8] k1_cn    chemical-potential CN-dependence mu = mu0*(1+k1cn*CN) (SI Eq. 84;
                  ES1-dominant, erf-fast decay, H/Cl signs opposite matching mu signs;
                  file value H 0.775 reproduces the measured stretch signal magnitude)
  L8[5] increment atomic core increment (moved the printed increments by exactly N_A*delta;
                  H's stored value is 0.0, matching Term 1's measurement)
  L8[6] kq2_rep   quadratic charge-dependence of Zeff (rep-only, q^2-scaled: big in NH4+)
  L9[*] ACP block (all slots move electronic/eps/gap only -- s,p,d,f projector pairs)
GLOBALS:
  G1[3] kpen1     penetration 1/R coefficient, elements WITH a core
  G1[8] kpen1_hhe penetration 1/R coefficient for H/He
  G2[0] kcn_glob  Eq.47 CN steepness (NEGATIVE -- the SI's Eq. 47 prints without the minus)
  G2[2] kpen2 ; G2[3] kpen3 ; G2[4] kpen4   (1/R^2..4 penetration terms)
  G1[9], G2[9]    dispersion globals (not yet named individually)
  G2[6], G2[7]    multipole electrostatics globals
Unnamed slots stay accessible by row/index; naming continues as terms get implemented.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "param", "gxtb", "gxtb_parameters")

import sys as _sys
_sys.path.insert(0, HERE)
import constants as _K
KEXP_REP = _K.REPULSION_KEXP    # from data/derived-constants.json (SI Sec. 1.6)


def parse(path=DEFAULT):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    glb, blocks, cur = [], {}, None
    for l in lines:
        t = l.split()
        if len(t) == 1 and re.fullmatch(r"\d+", t[0]):
            cur = int(t[0])
            blocks[cur] = []
        elif cur is None:
            glb.append([float(x) for x in t])
        else:
            blocks[cur].append([float(x) for x in t])
    g1, g2 = glb
    out = {"globals": {"g1": g1, "g2": g2,
                       "kpen1": g1[3], "kpen1_hhe": g1[8],
                       "kpen2": g2[2], "kpen3": g2[3], "kpen4": g2[4],
                       "kcn_glob": g2[0]},
           "element": {}}
    for z, rows in blocks.items():
        l1 = rows[0]
        out["element"][z] = {
            "alpha0": l1[0], "r0": l1[1], "zeff0": l1[2], "kcn_rep": l1[3],
            "kq_rep": l1[4], "rcov_cn": l1[5], "l1": l1,
            "shells": [rows[i] for i in range(1, 7)],       # L2..L7, 4-wide (s,p,d,f)
            "l8": rows[7], "l9": rows[8],
            "increment": rows[7][5], "kq2_rep": rows[7][6],
        }
    return out


if __name__ == "__main__":
    import json
    P = parse()
    print(f"globals: {len(P['globals']['g1'])}+{len(P['globals']['g2'])}; "
          f"elements: {len(P['element'])}")
    # cross-check the located increments against Term 1's oracle-measured values
    inc = json.load(open(os.path.join(HERE, "data", "increments-v1.json")))["increments_eh"]
    sym2z = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16,
             "Cl": 17, "Br": 35, "I": 53, "Pd": 46}
    worst = max(abs(P["element"][sym2z[s]]["increment"] - v) for s, v in inc.items())
    print(f"increments (file slot L8[5]) vs Term-1 oracle measurements, 13 elements: "
          f"max|diff| = {worst:.2e}  {'MATCH' if worst < 5e-7 else 'MISMATCH'}")
    assert worst < 5e-7
    print("PARAMS PARSER SELFTEST PASSED")
