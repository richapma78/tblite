"""Measure the ES3 target scalar: printed ES2+3 minus OUR ES2 (the validated engine
formulas -- onsite decoded gamma2 + offsite plain KO -- evaluated at the oracle density).
Controls first (F2, H2O: our ES2 matched to 1e-5/1e-3, so target ~ 0 there validates the
instrument); then the ions. Sign and size of the ion targets decide whether ES3 is the
missing ion layer or the hunt must redirect."""
import math
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import gxtb_engine as GE  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 2: "He", 8: "O", 9: "F"}


def tri(R):
    return [[0, 0, 0], [R, 0, 0], [R / 2, R * math.sqrt(3) / 2, 0]]


ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
CASES = [
    ("F2", [9, 9], [[0, 0, 0], [0, 0, 2.668]], 0),
    ("H2O", [8, 1, 1],
     [[0.0, 0.0, 0.0],
      [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
      [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]], 0),
    ("HeH+", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 1),
    ("H3+", [1, 1, 1], tri(1.65), 1),
    ("OH-", [8, 1], [[0, 0, 0], [0, 0, 1.83]], -1),
]

for name, zs, xyz, chg in CASES:
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=chg)
    P = rec["state"]["P"]
    B = GE.build(zs, np.array(xyz, float), charge=chg)
    S, meta, E, n = B["S"], B["meta"], B["E"], B["n"]
    m = np.diag((P / 2.0) @ S)
    q = {}
    for i in range(n):
        at, z, l, _ = meta[i]
        q[(at, l)] = q.get((at, l), 0.0) + 2 * m[i]
    for (at, l) in list(q):
        z = zs[at]
        q[(at, l)] = E[z]["ref"][l] - q[(at, l)]
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    es2 = 0.0
    for (at, l), qa in q.items():
        z = zs[at]
        for l2 in range(E[z]["nsh"]):
            g2 = GE.SRULE[z] * 2 * E[z]["U"][l] * E[z]["U"][l2] / \
                (E[z]["U"][l] + E[z]["U"][l2])
            es2 += 0.5 * qa * q[(at, l2)] * g2
    for (at, l), qa in q.items():
        for (bt, l2), qb in q.items():
            if bt <= at:
                continue
            z, zb = zs[at], zs[bt]
            gko = 1.0 / (float(B["Rab"][at, bt])
                         + 0.5 * (1.0 / E[z]["U"][l] + 1.0 / E[zb]["U"][l2]))
            es2 += qa * qb * gko
    p23 = float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)",
                          rec["state"]["raw"], re.M).group(1))
    qs = "  ".join(f"q({SYM[zs[a]]}{a})={qat[a]:+.3f}" for a in sorted(qat))
    print(f"{name:5s} printed ES2+3 {p23:+.5f}  ours-ES2 {es2:+.5f}  "
          f"TARGET(ES3?) {p23 - es2:+.5f}   [{qs}]")
