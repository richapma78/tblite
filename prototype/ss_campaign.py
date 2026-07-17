"""grand5: THE SS CAMPAIGN. The ss class is the one off-diagonal law that never closed
(rms 0.011, no cross-system transfer) and the HF LUMO sits on it. This round adds every
fully-invertible system that carries fresh ss physics: H3+ (pure ss, three-center), HeH+
and He2 (the banked He kernel), OH- (heteronuclear ss, new charge state). Refit ss on the
pooled data with per-family leave-one-out.

PRE-DECLARED INSTALL GATE (written before any run): the class law with the new ss law
installs ONLY if all SIX gate systems (H2, F2, HF, H2O, H3+, HeH+) hold worst printed-
window |eps d| <= 0.05 Eh. Anything less: banked, not installed, v1 stays. Per-system
attribution reported either way (the diagonal law's transfer to IONS is itself under
test here -- an ion failure may be the diag law's, not ss's; say which)."""
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
SYM = {1: "H", 2: "He", 8: "O", 9: "F"}


def tri(R):
    return [[0, 0, 0], [R, 0, 0], [R / 2, R * math.sqrt(3) / 2, 0]]


NEW = [
    ("H3+@1.5", [1, 1, 1], tri(1.5), 1, 2),
    ("H3+@1.65", [1, 1, 1], tri(1.65), 1, 2),
    ("H3+@1.9", [1, 1, 1], tri(1.9), 1, 2),
    ("HeH+@1.46", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 1, 2),
    ("HeH+@2.0", [2, 1], [[0, 0, 0], [0, 0, 2.0]], 1, 2),
    ("HeH+@2.5", [2, 1], [[0, 0, 0], [0, 0, 2.5]], 1, 2),
    ("He2@3.0", [2, 2], [[0, 0, 0], [0, 0, 3.0]], 0, 4),
    ("He2@4.0", [2, 2], [[0, 0, 0], [0, 0, 4.0]], 0, 4),
    ("OH-@1.83", [8, 1], [[0, 0, 0], [0, 0, 1.83]], -1, 10),
]


def extract(name, zs, xyz, charge):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge)
    P = rec["state"]["P"]
    B = GE.build(zs, np.array(xyz, float), charge=charge)
    GO, GC = GE.GRAND_O, GE.GRAND_C
    GE.GRAND_O, GE.GRAND_C = {}, None
    F_asm = GE.fock(P, B)
    GE.GRAND_O, GE.GRAND_C = GO, GC
    rem = rec["F"] - F_asm
    S, meta, n = B["S"], B["meta"], B["n"]
    mvec = np.diag((P / 2.0) @ S)
    gon = [2 * B["E"][meta[i][1]]["cx"] * B["E"][meta[i][1]]["L5"][meta[i][2]]
           for i in range(n)]
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            if meta[i][0] == meta[j][0]:
                continue
            cls = GE.pair_class(meta[i], meta[j], B["xyz"])
            out.append(dict(sys=name.split("@")[0], cls=cls, s=float(S[i, j]),
                            p=float(P[i, j]) / 2.0, gb=0.5 * (gon[i] + gon[j]),
                            gh=2 * gon[i] * gon[j] / (gon[i] + gon[j] + 1e-300),
                            mm=0.5 * float(mvec[i] + mvec[j]),
                            R=float(B["Rab"][meta[i][0], meta[j][0]]),
                            y=float(rem[i, j])))
    # diagonal remainders too: the diag law's ion transfer is reported (not refit)
    dmax = max(abs(float(rem[i, i])) for i in range(n))
    return out, dmax


D = json.load(open("/mnt/c/Projects/tblite-gxtb/prototype/data/grand-short.json"))
# pull the OLD ss rows out of the banked class dataset by re-deriving from offdiag +
# diag banks (they carry no gh) -- simpler: re-extract the four old systems too, same
# code path, so every row has identical provenance
OLD = []
for R in (1.2, 1.4, 1.7, 2.0, 2.5, 3.0):
    OLD.append((f"H2@{R}", [1, 1], [[0, 0, 0], [0, 0, R]], 0, 2))
for R in (2.668, 3.0, 3.4, 3.8):
    OLD.append((f"F2@{R}", [9, 9], [[0, 0, 0], [0, 0, R]], 0, 14))
for R in (1.733, 2.0, 2.4, 2.7):
    OLD.append((f"HF@{R}", [1, 9], [[0, 0, 0], [0, 0, R]], 0, 8))
ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
OLD.append(("H2O", [8, 1, 1],
            [[0.0, 0.0, 0.0],
             [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
             [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]], 0, 8))

rows = []
for name, zs, xyz, chg, _nel in OLD + NEW:
    try:
        rr, dmax = extract(name, zs, xyz, chg)
        rows.extend(rr)
        tag = "" if name.split("@")[0] in ("H2", "F2", "HF", "H2O") else \
            f"   [ION/He diag-law transfer: |rem_diag|max {dmax:.4f}]"
        print(f"  {name}: {len(rr)} cross rows{tag}", flush=True)
    except Exception as e:
        print(f"  {name}: FAILED ({type(e).__name__}: {e})", flush=True)

ss = [r for r in rows if r["cls"] == "ss"]
fams = sorted({r["sys"] for r in ss})
print(f"\nss rows: {len(ss)} across families {fams}")

FEATS = {
    "gb_s":   lambda r: r["gb"] * r["s"],
    "gb_s3":  lambda r: r["gb"] * r["s"] ** 3,
    "gh_s":   lambda r: r["gh"] * r["s"],
    "gh_s3":  lambda r: r["gh"] * r["s"] ** 3,
    "gb_p":   lambda r: r["gb"] * r["p"],
    "gb_p_s2": lambda r: r["gb"] * r["p"] * r["s"] ** 2,
    "s":      lambda r: r["s"],
    "p":      lambda r: r["p"],
    "obj":    lambda r: 0.42 * math.erf(0.2347 * r["R"]) / r["R"] - 0.048 * r["s"],
    "obj_p":  lambda r: (0.42 * math.erf(0.2347 * r["R"]) / r["R"]
                         - 0.048 * r["s"]) * r["p"],
    "gb_s_mm": lambda r: r["gb"] * r["s"] * r["mm"],
    "gh_p":   lambda r: r["gh"] * r["p"],
}
names = list(FEATS)
y = np.array([r["y"] for r in ss])
M = np.array([[FEATS[nm](r) for nm in names] for r in ss])
print(f"raw ss rms {float(np.sqrt((y**2).mean())):.4f}")
res = []
for k in (1, 2, 3):
    for combo in itertools.combinations(range(len(names)), k):
        A = M[:, list(combo)]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        if np.abs(c).max() > 50:
            continue
        rms = float(np.sqrt(np.mean((A @ c - y) ** 2)))
        res.append((rms, combo, c))
res.sort(key=lambda t: t[0])
print("top ss fits:")
for rms, combo, c in res[:5]:
    print(f"  rms {rms:.4f}: " + " ".join(f"{ci:+.4f}*{names[i]}"
                                          for ci, i in zip(c, combo)))
rms, combo, c = res[0]
print("\nper-family leave-one-out (best form):")
worst_loo = 0.0
for f in fams:
    tr = [k2 for k2, r in enumerate(ss) if r["sys"] != f]
    te = [k2 for k2, r in enumerate(ss) if r["sys"] == f]
    ct, *_ = np.linalg.lstsq(M[tr][:, list(combo)], y[tr], rcond=None)
    e = np.abs(M[te][:, list(combo)] @ ct - y[te])
    worst_loo = max(worst_loo, e.max())
    print(f"  hold out {f:5s}: |err|max {e.max():.4f}  rms {np.sqrt((e**2).mean()):.4f}")

D["ss_campaign"] = {"n": len(ss), "families": fams,
                    "best": {"terms": [names[i] for i in combo],
                             "coefs": [float(x) for x in c], "rms": rms},
                    "worst_loo": worst_loo}
json.dump(D, open("/mnt/c/Projects/tblite-gxtb/prototype/data/grand-short.json", "w"),
          indent=1)
print("\nbanked -> ss_campaign")
