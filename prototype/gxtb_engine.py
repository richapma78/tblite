"""gxtb_engine.py -- the GENERAL assembly engine: one code path for any closed-shell
molecule, every constant from the decoded laws. The Fortran port's blueprint.

Per-element inputs (all from the parameter file + decoded laws, NO per-molecule fits):
  levels -L2_l; mu_l; U_l; L5_l; c_x = s_rule(element)/9.59; metric k_l = L4_l (the
  L4-direct law; known bands: H +1.5%, F-s +0.7%, F-p +7%), kb = 0.4223*L8[3] (the kb law);
  kdiat sigma/pi = L8[0]/L8[1] (harmonic pair mean); kW_s/kW_p = G1[0]/G1[1]; ACP analytic.
Assembly: H0 (Eq-64 with arithmetic level means), ES1 (mu potentials), ES2 (onsite decoded
kernel + offsite KO, self-consistent Mulliken shell charges), X (per-AO same-shell
skeleton), unified Mulliken potential; short pieces/objects only where fitted (H2/F2/HF
pair tables) -- ABSENT for new pairs, so new-molecule gates are pure law-prediction.

First polyatomic gate: H2O (6 AOs, fully invertible). Pre-declared: all 6 eigenvalues
within 0.05 Eh (integration grade, no short pieces for O-H); report per orbital.
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
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

P_ = params.parse()
KW = [P_["globals"]["g1"][0], P_["globals"]["g1"][1]]
SRULE_P2 = (0.08247, 0.0920)
SRULE = {1: 0.4726, 2: 0.9218}
NVAL = {6: 4, 7: 5, 8: 6, 9: 7}
for z, nv in NVAL.items():
    SRULE[z] = SRULE_P2[0] * nv + SRULE_P2[1]
Bq = basisq.parse()


def elem(z):
    e = P_["element"][z]
    nsh = len(Bq[z]["shells"])
    return {"L2": e["shells"][0][:nsh], "MU": e["shells"][5][:nsh],
            "U": e["shells"][4][:nsh], "L5": e["shells"][3][:nsh],
            "cx": SRULE[z] / 9.59, "k": [e["shells"][2][l] for l in range(nsh)],
            "kb": 0.4223 * e["l8"][3], "kd_sg": e["l8"][0], "kd_pi": e["l8"][1],
            "ref": K.REFOCC[z], "nsh": nsh}


def build(zs, xyz_bohr, charge=0):
    E = {z: elem(z) for z in set(zs)}
    xyz = np.array(xyz_bohr, float)
    S = overlap.overlap(list(zs), xyz, charge=charge, ao_order="oracle")
    cns = adapt.basis_cn(list(zs), xyz)
    sh, meta = [], []
    for at, z in enumerate(zs):
        q_ = E[z]["kb"] * math.sqrt(cns[at])
        for l, prims in Bq[z]["shells"]:
            e_ = np.array([p[0] for p in prims])
            c0 = np.array([p[1] for p in prims])
            c1 = np.array([p[2] for p in prims])
            sh.append({"l": l, "at": at, "exp": e_ * E[z]["k"][l],
                       "coef": (c0 + c1 * q_).copy()})
            for m in range(2 * l + 1):
                meta.append((at, z, l, m))
    Sh = overlap.overlap(list(zs), xyz, shells=sh, ao_order="oracle")
    n = len(meta)
    H0 = np.zeros((n, n))
    for i in range(n):
        H0[i, i] = -E[meta[i][1]]["L2"][meta[i][2]]
    for i in range(n):
        for j in range(n):
            ai, zi, li, _mi = meta[i]
            aj, zj, lj, _mj = meta[j]
            if ai == aj:
                continue
            # sigma/pi resolution: project the AO pair onto the bond axis is already in
            # S-tilde; kdiat by the dominant type -- s-involving pairs are sigma; p-p pairs
            # mix sigma and pi per orientation, carried by S~'s frame content. v1: use the
            # HARMONIC sigma constant for s-any, and for p-p the sigma constant too --
            # EXCEPT pure-pi geometry (handled adequately for planar/diatomic frames by
            # the overlap's own vanishing sigma content). Known v1 limitation, labeled.
            kdi = 2 * E[zi]["kd_sg"] * E[zj]["kd_sg"] / (E[zi]["kd_sg"] + E[zj]["kd_sg"])
            if li == 1 and lj == 1:
                kpi = 2 * E[zi]["kd_pi"] * E[zj]["kd_pi"] / (E[zi]["kd_pi"] + E[zj]["kd_pi"])
                kdi = kpi if abs(Sh[i, j]) < 0.02 else kdi
            a = (KW[li] + KW[lj]) / 2 * kdi
            h = -(E[zi]["L2"][li] + E[zj]["L2"][lj]) / 2
            H0[i, j] = a * h * float(Sh[i, j])
    A = f2_stretch.acp_matrix(list(zs), xyz, charge=charge)
    nat = len(zs)
    Rab = np.zeros((nat, nat))
    for a in range(nat):
        for b in range(nat):
            Rab[a, b] = float(np.linalg.norm(xyz[a] - xyz[b]))
    return {"S": S, "H0": H0, "A": A, "meta": meta, "E": E, "zs": zs, "n": n, "Rab": Rab}


def fock(P, B):
    S, meta, E = B["S"], B["meta"], B["E"]
    n = B["n"]
    Ps = P / 2.0
    m = np.diag(Ps @ S)
    q = {}
    for i in range(n):
        at, z, l, _ = meta[i]
        q[(at, l)] = q.get((at, l), 0.0) + 2 * m[i]
    for (at, l) in list(q):
        z = B["zs"][at]
        q[(at, l)] = E[z]["ref"][l] - q[(at, l)]
    v = np.zeros(n)
    Rab = B["Rab"]
    for i in range(n):
        at, z, l, _ = meta[i]
        v[i] += E[z]["MU"][l]
        for l2 in range(E[z]["nsh"]):
            g2 = SRULE[z] * 2 * E[z]["U"][l] * E[z]["U"][l2] / \
                (E[z]["U"][l] + E[z]["U"][l2])
            v[i] += g2 * q[(at, l2)]
        # offsite ES2: the plain Klopman-Ohno kernel (k2x = 0; validated out-of-sample
        # on HF), generalized to the N-atom sum
        for bt in range(len(B["zs"])):
            if bt == at:
                continue
            zb = B["zs"][bt]
            for l2 in range(E[zb]["nsh"]):
                gko = 1.0 / (Rab[at, bt] + 0.5 * (1.0 / E[z]["U"][l]
                                                  + 1.0 / E[zb]["U"][l2]))
                v[i] += gko * q[(bt, l2)]
        v[i] += 2 * E[z]["cx"] * E[z]["L5"][l] * m[i]
    return B["H0"] + B["A"] - 0.5 * S * (v[:, None] + v[None, :])


def scf(zs, xyz_bohr, nel, charge=0, iters=120, mix=0.4):
    B = build(zs, xyz_bohr, charge)
    S = B["S"]
    s, Uo = np.linalg.eigh(S)
    X = Uo @ np.diag(1 / np.sqrt(s)) @ Uo.T
    nocc = nel // 2
    P = np.zeros((B["n"], B["n"]))
    for _ in range(iters):
        F = fock(P, B)
        w, Vp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Vp
        Pn = 2.0 * C[:, :nocc] @ C[:, :nocc].T
        if np.abs(Pn - P).max() < 1e-9:
            P = Pn
            break
        P = mix * Pn + (1 - mix) * P
    F = fock(P, B)
    w, _ = np.linalg.eigh(X.T @ F @ X)
    return w, P


def main():
    # H2O: the first polyatomic, PURE LAW-PREDICTION (no fitted short pieces for O-H)
    ang = 104.5 * math.pi / 180
    r_oh = 0.9572 * K.BOHR
    xyz = [[0.0, 0.0, 0.0],
           [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
           [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
    zs = [8, 1, 1]
    w, P = scf(zs, xyz, nel=8)
    atoms = [("O", 0, 0, 0),
             ("H", xyz[1][0] / K.BOHR, 0, xyz[1][2] / K.BOHR),
             ("H", xyz[2][0] / K.BOHR, 0, xyz[2][2] / K.BOHR)]
    r = oracle.run(atoms)
    ref = sorted(x / K.EV for x in r["eps_ev"])
    d = [a - b for a, b in zip(sorted(w), ref)]
    worst = max(abs(x) for x in d)
    B = build(zs, np.array(xyz))
    m = np.diag((P / 2.0) @ B["S"])
    qO = sum(K.REFOCC[8].values()) - 2 * sum(m[i] for i in range(B["n"])
                                             if B["meta"][i][0] == 0)
    print(f"H2O v2 (+offsite-ES2): worst |d| = {worst:.4f}  "
          f"{'PASS' if worst <= 0.05 else 'MISS'} (bar 0.05; O-H short pieces + mu-CN "
          f"still absent)   q(O) = {qO:+.3f}")
    for a, b, dd in zip(sorted(w), ref, d):
        print(f"    {a:+.4f}  vs {b:+.4f}   d {dd:+.4f}")


if __name__ == "__main__":
    main()
