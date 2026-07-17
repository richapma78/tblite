"""ES3 install verification.
1. Solve OH-'s tau entries (1-8/8-1 @ 1.83) from its pass-76 Euler cleans and write
   them into the banked table.
2. ION REMAINDER RE-EXTRACTION (pre-declared): with ES3 installed, the diagonal Fock
   remainders of HeH+@1.46, H3+@1.65, OH-@1.83 must drop below 0.05 (from 0.212 /
   0.042 / 0.194). This tests the POTENTIAL (v3 = dE3/dq Mulliken shift) -- the
   energies are exact-by-construction at these geometries.
3. SIX-SYSTEM GATE: H2, F2, HF, H2O, H3+, HeH+ full SCF vs oracle eigenvalues,
   bar 0.05 per system. (H2O/HF carry a labeled double-count caveat: the grand diag
   law absorbed some ES3; the re-fit on clean remainders is the next push.)"""
import importlib
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import oracle  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 2: "He", 8: "O", 9: "F"}

# ---- 1. OH- tau solve (Euler cleans from pass 76: G_H -0.00527, G_O +0.01258)
zs = [8, 1]
xyz = [[0, 0, 0], [0, 0, 1.83]]
atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
         for z, (x, y_, zc) in zip(zs, xyz)]
rec = fock_recon.fock_ao(atoms, charge=-1)
P, S = rec["state"]["P"], rec["S"]
m = np.diag((P / 2.0) @ S)
meta = [(0, 0), (0, 1), (0, 1), (0, 1), (1, 0)]
q = {}
for mm, (at, l) in zip(m, meta):
    q[(at, l)] = q.get((at, l), 0.0) + 2 * mm
for (at, l) in list(q):
    q[(at, l)] = K.REFOCC[zs[at]].get(l, 0.0) - q[(at, l)]
qat = {}
for (at, l), v in q.items():
    qat[at] = qat.get(at, 0.0) + v
GAM = {8: 0.1454, 1: 0.8141717488}
c_O = sum((1 / 6) * qa * qb * qat[at] * GAM[8]
          for (at, la), qa in q.items() if at == 0
          for (bt, lb), qb in q.items() if bt == 1)
c_H = sum((1 / 6) * qa * qb * qat[at] * GAM[1]
          for (at, la), qa in q.items() if at == 1
          for (bt, lb), qb in q.items() if bt == 0)
tau_O = +0.01258 / c_O
tau_H = -0.00527 / c_H
print(f"OH- taus @1.83: tau_O = {tau_O:+.5f}  tau_H = {tau_H:+.5f}")
p = "/mnt/c/Projects/tblite-gxtb/prototype/data/derived-constants.json"
d = json.load(open(p))
d["es3"]["tau_pairs"]["8-1"] = [[1.809, -0.12710], [1.83, round(tau_O, 5)]]
d["es3"]["tau_pairs"]["1-8"] = [[1.809, 0.20980], [1.83, round(tau_H, 5)]]
json.dump(d, open(p, "w"), indent=1)
importlib.reload(GE)
print("banked + engine reloaded")

# ---- 2. ion diagonal remainders with ES3 installed
print("\nION REMAINDER RE-EXTRACTION (pre-declared bar 0.05):")


def tri(R):
    return [[0, 0, 0], [R, 0, 0], [R / 2, R * math.sqrt(3) / 2, 0]]


for name, zs_, xyz_, chg, prev in (
        ("HeH+@1.46", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 1, 0.212),
        ("H3+@1.65", [1, 1, 1], tri(1.65), 1, 0.042),
        ("OH-@1.83", [8, 1], [[0, 0, 0], [0, 0, 1.83]], -1, 0.194)):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs_, xyz_)]
    r2 = fock_recon.fock_ao(atoms, charge=chg)
    B = GE.build(zs_, np.array(xyz_, float), charge=chg)
    F_asm = GE.fock(r2["state"]["P"], B)
    rem = r2["F"] - F_asm
    dmax = max(abs(float(rem[i, i])) for i in range(B["n"]))
    print(f"  {name}: |rem_diag|max {dmax:.4f} (was {prev})  "
          f"{'PASS' if dmax < 0.05 else 'MISS'}", flush=True)

# ---- 3. six-system gate
print("\nSIX-SYSTEM GATE (bar 0.05, printed window):")
ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
wxyz = [[0.0, 0.0, 0.0],
        [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
        [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
cases = [("H2", [1, 1], [[0, 0, 0], [0, 0, 1.4]], 2, 0),
         ("F2", [9, 9], [[0, 0, 0], [0, 0, 2.668]], 14, 0),
         ("HF", [1, 9], [[0, 0, 0], [0, 0, 1.733]], 8, 0),
         ("H2O", [8, 1, 1], wxyz, 8, 0),
         ("H3+", [1, 1, 1], tri(1.65), 2, 1),
         ("HeH+", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 2, 1)]
npass = 0
for name, zs_, xyz_, nel, chg in cases:
    w, P = GE.scf(zs_, xyz_, nel=nel, charge=chg)
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs_, xyz_)]
    r = oracle.run(atoms, charge=chg)
    ref = sorted(x / K.EV for x in r["eps_ev"])
    dd = [a - b for a, b in zip(sorted(w)[:len(ref)], ref)]
    worst = max(abs(x) for x in dd)
    ok = worst <= 0.05
    npass += ok
    print(f"  {name:5s} worst |d| = {worst:.4f}  {'PASS' if ok else 'MISS'}   "
          + " ".join(f"{x:+.3f}" for x in dd), flush=True)
print(f"\nSIX-SYSTEM GATE: {npass}/6")
