"""grand4: the SIGMA/PI class-resolved off-diagonal round. Re-extract the cross-atom
remainders with class labels (gxtb_engine.pair_class -- ONE classifier shared with the
engine), fit each class independently (small feature family; OBJ candidates allowed where
the class is sigma-involving), parsimony rule: take the 1-term fit unless 2 terms improve
class-rms by >15%. Bank as offdiag_class. PRE-DECLARED INSTALL GATE (recorded before the
run): installed, the four-system gate must stay >= 3/4 with no system's worst |d| over
0.05; SUCCESS = H2O LUMO under the bar (4/4); any regression -> banked, NOT installed."""
import itertools
import json
import math
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import gxtb_engine as GE  # noqa: E402

BOHR = K.BOHR


def systems():
    out = []
    for R in (1.2, 1.4, 1.7, 2.0, 2.5, 3.0):
        out.append((f"H2@{R}", [1, 1], [[0, 0, 0], [0, 0, R]]))
    for R in (2.668, 3.0, 3.4, 3.8):
        out.append((f"F2@{R}", [9, 9], [[0, 0, 0], [0, 0, R]]))
    for R in (1.733, 2.0, 2.4, 2.7):
        out.append((f"HF@{R}", [1, 9], [[0, 0, 0], [0, 0, R]]))
    ang = 104.5 * math.pi / 180
    r_oh = 0.9572 * BOHR
    out.append(("H2O", [8, 1, 1],
                [[0.0, 0.0, 0.0],
                 [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
                 [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]))
    return out


SYM = {1: "H", 8: "O", 9: "F"}
rows = []
for name, zs, xyz in systems():
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms)
    P = rec["state"]["P"]
    B = GE.build(zs, np.array(xyz, float))
    # remainder against the engine WITHOUT its old flat off-diag law: strip it so the
    # class fit is not fitting on top of the very layer it replaces
    GO_save, GC_save = GE.GRAND_O, GE.GRAND_C
    GE.GRAND_O, GE.GRAND_C = {}, None
    F_asm = GE.fock(P, B)
    GE.GRAND_O, GE.GRAND_C = GO_save, GC_save
    rem = rec["F"] - F_asm
    S, meta, n = B["S"], B["meta"], B["n"]
    mvec = np.diag((P / 2.0) @ S)
    gon = [2 * B["E"][meta[i][1]]["cx"] * B["E"][meta[i][1]]["L5"][meta[i][2]]
           for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if meta[i][0] == meta[j][0]:
                continue
            cls = GE.pair_class(meta[i], meta[j], B["xyz"])
            R = float(B["Rab"][meta[i][0], meta[j][0]])
            rows.append(dict(sys=name.split("@")[0], cls=cls, s=float(S[i, j]),
                             p=float(P[i, j]) / 2.0, gb=0.5 * (gon[i] + gon[j]),
                             mm=0.5 * float(mvec[i] + mvec[j]), R=R, y=float(rem[i, j])))

FEATS = {
    "gb_s":   lambda r: r["gb"] * r["s"],
    "gb_s3":  lambda r: r["gb"] * r["s"] ** 3,
    "gb_p":   lambda r: r["gb"] * r["p"],
    "gb_p_s2": lambda r: r["gb"] * r["p"] * r["s"] ** 2,
    "s":      lambda r: r["s"],
    "p":      lambda r: r["p"],
    "obj":    lambda r: 0.42 * math.erf(0.2347 * r["R"]) / r["R"] - 0.048 * r["s"],
    "obj_p":  lambda r: (0.42 * math.erf(0.2347 * r["R"]) / r["R"]
                         - 0.048 * r["s"]) * r["p"],
    "gb_s_mm": lambda r: r["gb"] * r["s"] * r["mm"],
}
names = list(FEATS)
by_cls = {}
for r in rows:
    by_cls.setdefault(r["cls"], []).append(r)
print("class populations: " + "  ".join(f"{c}:{len(v)}" for c, v in sorted(by_cls.items())))

laws = {}
for cls, rs in sorted(by_cls.items()):
    y = np.array([r["y"] for r in rs])
    raw = float(np.sqrt((y**2).mean()))
    if raw < 2e-3:
        laws[cls] = {"terms": [], "coefs": [], "rms": raw, "n": len(rs)}
        print(f"  {cls:5s} n={len(rs):3d} raw rms {raw:.4f} -> ZERO law (below noise)")
        continue
    M = np.array([[FEATS[nm](r) for nm in names] for r in rs])
    # STATIC restriction where declared: density-carried off-diag terms create SCF
    # feedback (measured: HF LUMO overshoot at self-consistency with a remainder fit
    # that is exact at the oracle density). sp_s is restricted to static forms.
    allowed = names
    allowed_idx = [names.index(a) for a in allowed]
    best = {}
    kmax = 3 if len(rs) >= 12 else 2
    for k in range(1, kmax + 1):
        cands = []
        for combo in itertools.combinations(allowed_idx, k):
            A = M[:, list(combo)]
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            if np.abs(c).max() > 50:      # collinearity guard
                continue
            rms = float(np.sqrt(np.mean((A @ c - y) ** 2)))
            cands.append((rms, combo, c))
        if cands:
            best[k] = min(cands, key=lambda t: t[0])
    rms, combo, c = best[1]
    for k in range(2, kmax + 1):
        if k in best and len(rs) >= 3 * k and best[k][0] < 0.85 * rms:
            rms, combo, c = best[k]
    laws[cls] = {"terms": [names[i] for i in combo], "coefs": [float(x) for x in c],
                 "rms": rms, "n": len(rs)}
    print(f"  {cls:5s} n={len(rs):3d} raw rms {raw:.4f} -> {rms:.4f}  "
          + " ".join(f"{ci:+.4f}*{names[i]}" for ci, i in zip(c, combo)))

# H2 leave-out check on the ss class (H2 is all-ss; can the others predict it?)
ss = by_cls.get("ss", [])
tr = [r for r in ss if r["sys"] != "H2"]
te = [r for r in ss if r["sys"] == "H2"]
if tr and te and laws["ss"]["terms"]:
    combo = [names.index(t) for t in laws["ss"]["terms"]]
    Mtr = np.array([[FEATS[names[i]](r) for i in combo] for r in tr])
    ct, *_ = np.linalg.lstsq(Mtr, np.array([r["y"] for r in tr]), rcond=None)
    Mte = np.array([[FEATS[names[i]](r) for i in combo] for r in te])
    e = np.abs(Mte @ ct - np.array([r["y"] for r in te]))
    print(f"\nss-class H2-holdout: |err|max {e.max():.4f}  rms {np.sqrt((e**2).mean()):.4f}")

D = json.load(open("/mnt/c/Projects/tblite-gxtb/prototype/data/grand-short.json"))
D["offdiag_class"] = {"classes": laws, "installed": False,
                      "note": "sigma/pi class-resolved; ZERO law = class below noise. "
                              "NOT INSTALLED (pre-declared gate): dominates v1 on "
                              "H2/F2/H2O (water closes fully, 0.017) but the HF LUMO "
                              "regresses 0.026->0.053, 0.003 over the bar -- it sits on "
                              "the unresolved ss channel (rms 0.011, H2-holdout fails). "
                              "Static-restriction variants of sp_s and ss both REFUTED "
                              "(HF worse: 0.056/0.057). Next entry ticket: crack ss."}
json.dump(D, open("/mnt/c/Projects/tblite-gxtb/prototype/data/grand-short.json", "w"),
          indent=1)
print("banked -> offdiag_class")
