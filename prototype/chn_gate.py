"""THE PORT-SCOPE GATE: carbon and nitrogen, never before tested. CH4 and NH3 full
SCF vs oracle eigenvalues, bar 0.05. Every law extrapolates from its H/He/O/F fitting
domain via the element-generic structure (constants from the parameter file). The
result -- pass or miss -- defines the Fortran port's validated scope; no tuning
either way."""
import math
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import oracle  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N"}

a = 1.087 * BOHR / math.sqrt(3)
CH4 = ([6, 1, 1, 1, 1],
       [[0, 0, 0], [a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]], 8)
r_nh = 1.012 * BOHR
st, ct = 0.9262, -0.3770
NH3 = ([7, 1, 1, 1],
       [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                       r_nh * st * math.sin(2 * math.pi * k / 3),
                       r_nh * ct] for k in range(3)], 8)

for name, (zs, xyz, nel) in (("CH4", CH4), ("NH3", NH3)):
    w, P = GE.scf(zs, xyz, nel=nel)
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    r = oracle.run(atoms)
    ref = sorted(x / K.EV for x in r["eps_ev"])
    dd = [x - y for x, y in zip(sorted(w)[:len(ref)], ref)]
    worst = max(abs(x) for x in dd)
    occ = nel // 2
    worst_occ = max(abs(x) for x in dd[:occ])
    print(f"{name}: worst |d| = {worst:.4f} ({'PASS' if worst <= 0.05 else 'MISS'})  "
          f"occupied-only {worst_occ:.4f}"
          f" ({'PASS' if worst_occ <= 0.05 else 'MISS'})")
    print("   " + " ".join(f"{x:+.3f}" for x in dd))
