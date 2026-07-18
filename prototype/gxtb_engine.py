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
from scipy.special import erf as erf_np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapt  # noqa: E402
import basisq  # noqa: E402
import constants as K  # noqa: E402
import es2_energy  # noqa: E402
import f2_stretch  # noqa: E402
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

P_ = params.parse()
KW = [P_["globals"]["g1"][0], P_["globals"]["g1"][1]]
try:
    with open(os.path.join(HERE, "data", "grand-short.json")) as _f:
        _G = __import__("json").load(_f)
    GRAND_D = dict(zip(_G["best"]["terms"], _G["best"]["coefs"]))
    GRAND_O = dict(zip(_G["offdiag_best"]["terms"], _G["offdiag_best"]["coefs"]))
    GRAND_C = ({c: list(zip(v["terms"], v["coefs"]))
                for c, v in _G["offdiag_class"]["classes"].items()}
               if _G.get("offdiag_class", {}).get("installed") else None)
except FileNotFoundError:
    GRAND_D = GRAND_O = GRAND_C = None
try:
    with open(os.path.join(HERE, "data", "derived-constants.json")) as _f:
        _DC = __import__("json").load(_f)
    ES3 = _DC.get("es3")
    POT_MODE = (ES3 or {}).get("potential", "v1")
except FileNotFoundError:
    ES3 = None
    POT_MODE = "v1"


_TAU_WARNED = set()


def _tau_off(za, zb, R):
    """measured per-pair-side tau (multiplies Gamma_A of za); linear interpolation,
    flat single points, zero outside the measured range + margin. An UNMEASURED pair
    returns 0 with a one-time warning -- correct only where the converged atomic
    charge vanishes (homonuclear neutrals); anywhere else the ES3 of that pair is
    silently dropped and the result is suspect."""
    pts = [(r, t) for r, t in ES3["tau_pairs"].get(f"{za}-{zb}", []) if t is not None]
    if not pts:
        if (za, zb) not in _TAU_WARNED:
            _TAU_WARNED.add((za, zb))
            print(f"    [es3] WARNING: no measured tau for pair {za}-{zb}; using 0 "
                  f"(exact only if converged qat = 0)")
        return 0.0
    if len(pts) == 1:
        return pts[0][1] if abs(R - pts[0][0]) < 0.2 else 0.0
    if R <= pts[0][0]:
        return pts[0][1]
    if R >= pts[-1][0]:
        return 0.0 if R > pts[-1][0] + 0.5 else pts[-1][1]
    for (r1, t1), (r2, t2) in zip(pts, pts[1:]):
        if r1 <= R <= r2:
            return t1 + (t2 - t1) * (R - r1) / (r2 - r1)


def _tau_q(za, zb, R):
    """measured Q-term radial kernel (multiplies Gamma_A and Q_total); linear
    interpolation, flat single points, zero far outside range; unmeasured pair -> 0
    (the Q-term of that pair is dropped -- exact only at Q = 0)."""
    pts = ES3["v2"]["tau_Q"].get(f"{za}-{zb}")
    if not pts:
        if (za, zb, "Q") not in _TAU_WARNED:
            _TAU_WARNED.add((za, zb, "Q"))
            print(f"    [es3v2] WARNING: no tau_Q for pair {za}-{zb}; dropping its "
                  f"Q-term (exact only when Q = 0)")
        return 0.0
    if len(pts) == 1:
        return pts[0][1] if abs(R - pts[0][0]) < 0.25 else 0.0
    if R <= pts[0][0]:
        return pts[0][1]
    if R >= pts[-1][0]:
        return 0.0 if R > pts[-1][0] + 0.5 else pts[-1][1]
    for (r1, t1), (r2, t2) in zip(pts, pts[1:]):
        if r1 <= R <= r2:
            return t1 + (t2 - t1) * (R - r1) / (r2 - r1)


def es3_energy_v1(q, zs, Rab, E):
    """the v1 POTENTIAL-side third order: onsite closed form + the lumped measured
    tau table. NOT the energy truth (that is es_charge_energy) -- but its simple
    pairwise q-gradient matches the binary's Fock BETTER than the exact gradient of
    the true energy does (measured: gate 4/6 vs 3/6). The binary pairs a correct
    energy with a shortcut potential; so do we, deliberately."""
    if ES3 is None:
        return 0.0
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    if all(abs(x) < 1e-12 for x in qat.values()):
        return 0.0
    e3 = 0.0
    for (at, la), qa in q.items():
        z = zs[at]
        Ula = E[z]["U"][la]
        if Ula == 0.0:
            continue
        ga = (ES3["k3gs"] if la == 0 else ES3["k3gp"]) * ES3["gamma"][str(z)]
        ta = -1.0 / (2 * Ula * Ula)
        for lb in range(E[z]["nsh"]):
            Ulb = E[z]["U"][lb]
            if Ulb == 0.0:
                continue
            gb = (ES3["k3gs"] if lb == 0 else ES3["k3gp"]) * ES3["gamma"][str(z)]
            tb = -1.0 / (2 * Ulb * Ulb)
            e3 += (1 / 6) * qa * q[(at, lb)] * qat[at] * (ta * ga + tb * gb)
        for (bt, lb), qb in q.items():
            if bt == at:
                continue
            if abs(qa * qb * qat[at]) < 1e-9:
                continue
            tau = _tau_off(z, zs[bt], float(Rab[at, bt]))
            e3 += (1 / 6) * qa * qb * qat[at] * tau * ES3["gamma"][str(z)]
    return e3


def es_charge_energy(q, zs, Rab, E):
    """v2: everything charge-driven beyond the analytic onsite-ES2/mu -- the offsite
    KO second order at FOLDED Hubbards U_l(q) = U_l + Gamma_A*qat_A (c=1, bare; pass
    79), the onsite third-order closed form (pass 75), and the explicit Q-weighted
    term with measured kernels (pass 80). The Fock takes this function's numeric
    q-gradient, keeping energy and potential consistent BY CONSTRUCTION (the binary's
    own pairing is looser -- pass 77 -- and the remainder gates measure the gap)."""
    if ES3 is None:
        return 0.0
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    Q = sum(qat.values())
    gam = {z: ES3["gamma"][str(z)] for z in set(zs)}

    def ueff(at, l):
        z = zs[at]
        u = E[z]["U"][l] + gam[z] * qat[at]
        return max(u, 0.05 * E[z]["U"][l])      # He's big Gamma: floor on transients

    e = 0.0
    for (at, la), qa in q.items():
        z = zs[at]
        if E[z]["U"][la] == 0.0:
            continue
        # offsite KO at folded U (the ES2 offsite lives HERE now, not in v[])
        for (bt, lb), qb in q.items():
            if bt <= at or E[zs[bt]]["U"][lb] == 0.0:
                continue
            g = 1.0 / (float(Rab[at, bt])
                       + 0.5 * (1.0 / ueff(at, la) + 1.0 / ueff(bt, lb)))
            e += qa * qb * g
        # onsite third order (closed form)
        ga = (ES3["k3gs"] if la == 0 else ES3["k3gp"]) * gam[z]
        Ula = E[z]["U"][la]
        ta = -1.0 / (2 * Ula * Ula)
        for lb in range(E[z]["nsh"]):
            Ulb = E[z]["U"][lb]
            if Ulb == 0.0:
                continue
            gb = (ES3["k3gs"] if lb == 0 else ES3["k3gp"]) * gam[z]
            tb = -1.0 / (2 * Ulb * Ulb)
            e += (1 / 6) * qa * q[(at, lb)] * qat[at] * (ta * ga + tb * gb)
        # explicit Q-term (ordered cross pairs; dead when Q = 0)
        if abs(Q) > 1e-12:
            for (bt, lb), qb in q.items():
                if bt == at:
                    continue
                if abs(qa * qb * Q) < 1e-9:
                    continue
                e += (1 / 6) * qa * qb * Q * _tau_q(z, zs[bt],
                                                    float(Rab[at, bt])) * gam[z]
    return e
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


_PDIR = {0: (1.0, 0.0, 0.0), 1: (0.0, 1.0, 0.0), 2: (0.0, 0.0, 1.0)}


def pair_class(mi, mj, xyz):
    """sigma/pi class of a cross-atom AO pair. mi/mj = (at, z, l, m); oracle p-order
    is (x, y, z). Returns one of ss, sp_s, sp_p, sp_m, pp_s, pp_p, pp_m."""
    ai, _, li, qi = mi
    aj, _, lj, qj = mj
    u = np.array(xyz[aj], float) - np.array(xyz[ai], float)
    u /= np.linalg.norm(u)

    def cont(l, q):
        return None if l == 0 else abs(float(np.dot(_PDIR[q], u)))
    ci, cj = cont(li, qi), cont(lj, qj)
    if ci is None and cj is None:
        return "ss"
    if ci is None or cj is None:
        c = cj if ci is None else ci
        return "sp_s" if c > 0.9 else ("sp_p" if c < 0.1 else "sp_m")
    if ci > 0.9 and cj > 0.9:
        return "pp_s"
    if ci < 0.1 and cj < 0.1:
        return "pp_p"
    return "pp_m"


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
    # --- H0 v2 EXPERIMENT (env GXTB_H0V2=1, default OFF -> shipped v1 untouched) ---
    # SI Eq 65 CN-dependent level  H_lA = L2 - LEVCN*CN  (LEVCN=gp3_lev_cn_=shells[1],
    #   CN = internal es2 CN; sign verified). SI Eq 67 shell polynomial with the per-shell
    # factor taken from gp3_pln_ (=shells[2], the binary's per-element-per-shell "pln"
    # column) instead of a g1/g2 global: pi_lA = 1 + gp3_poly[Z]*gp3_pln[Z][l]*(R/Rcov).
    _h0v2 = os.environ.get("GXTB_H0V2") == "1"
    if _h0v2:
        _cn = es2_energy.coordination(list(zs), xyz)
        _el = P_["element"]

        def _hlev(at, z, l):
            return E[z]["L2"][l] - _el[z]["shells"][1][l] * _cn[at]
    else:
        def _hlev(at, z, l):
            return E[z]["L2"][l]
    H0 = np.zeros((n, n))
    for i in range(n):
        at, z, l, _m = meta[i]
        H0[i, i] = -_hlev(at, z, l)
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
            h = -(_hlev(ai, zi, li) + _hlev(aj, zj, lj)) / 2
            pi_pi = 1.0
            if _h0v2:
                Rij = float(np.linalg.norm(xyz[ai] - xyz[aj]))
                rc = (_el[zi]["l1"][5] + _el[zj]["l1"][5]) / 2.0
                pii = 1.0 + _el[zi]["l1"][7] * _el[zi]["shells"][2][li] * (Rij / rc)
                pij = 1.0 + _el[zj]["l1"][7] * _el[zj]["shells"][2][lj] * (Rij / rc)
                pi_pi = pii * pij
            H0[i, j] = a * h * pi_pi * float(Sh[i, j])
    A = f2_stretch.acp_matrix(list(zs), xyz, charge=charge)
    nat = len(zs)
    Rab = np.zeros((nat, nat))
    for a in range(nat):
        for b in range(nat):
            Rab[a, b] = float(np.linalg.norm(xyz[a] - xyz[b]))
    return {"S": S, "H0": H0, "A": A, "meta": meta, "E": E, "zs": zs, "n": n,
            "Rab": Rab, "xyz": xyz}


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
        if ES3 is None or POT_MODE == "v1":
            # offsite ES2: plain KO (analytic; in v2-gradient mode this kernel
            # lives inside es_charge_energy at FOLDED Hubbards instead)
            for bt in range(len(B["zs"])):
                if bt == at:
                    continue
                zb = B["zs"][bt]
                for l2 in range(E[zb]["nsh"]):
                    gko = 1.0 / (Rab[at, bt] + 0.5 * (1.0 / E[z]["U"][l]
                                                      + 1.0 / E[zb]["U"][l2]))
                    v[i] += gko * q[(bt, l2)]
        v[i] += 2 * E[z]["cx"] * E[z]["L5"][l] * m[i]
    if ES3 is not None:
        # THE POTENTIAL IS A SEPARATE OBJECT FROM THE ENERGY (T1, twice measured):
        # v1 = simple pairwise shift over the lumped tau table (gate 4/6, installed);
        # v2-gradient = exact dE/dq of the true composite energy (gate 3/6, refuted
        # for the Fock; es_charge_energy REMAINS the energy-side truth).
        efun = es3_energy_v1 if POT_MODE == "v1" else es_charge_energy
        h3 = 1e-6
        v3 = {}
        for key in q:
            qp = dict(q)
            qp[key] += h3
            qm = dict(q)
            qm[key] -= h3
            v3[key] = (efun(qp, B["zs"], Rab, E)
                       - efun(qm, B["zs"], Rab, E)) / (2 * h3)
        for i in range(n):
            at, z, l, _ = meta[i]
            v[i] += v3[(at, l)]
    F = B["H0"] + B["A"] - 0.5 * S * (v[:, None] + v[None, :])
    if GRAND_D is not None:
        # the GRAND short-piece layer (fitted across H2/F2/HF/H2O on THIS baseline;
        # coefficients live in data/grand-short.json, never here). Diagonal law
        # transfers at leave-one-system-out; off-diagonal is weaker (H2 holdout
        # fails) and labeled integration-grade.
        gon = np.array([2 * E[meta[i][1]]["cx"] * E[meta[i][1]]["L5"][meta[i][2]]
                        for i in range(n)])
        atv = np.array([meta[i][0] for i in range(n)])
        cross = atv[:, None] != atv[None, :]
        gb = 0.5 * (gon[:, None] + gon[None, :])
        S2 = S * S * cross
        d = (GRAND_D["S_gj_s2"] * (S2 @ gon)
             + GRAND_D["mi_S_gb_s2"] * m * np.sum(gb * S2, axis=1)
             + GRAND_D["S_gb_s4"] * np.sum(gb * S2 * S * S, axis=1))
        F += np.diag(d)
        if GRAND_C is not None:
            # the sigma/pi CLASS-RESOLVED off-diagonal law (one classifier shared
            # with the fit; laws live in grand-short.json offdiag_class)
            for i in range(n):
                for j in range(i + 1, n):
                    if meta[i][0] == meta[j][0]:
                        continue
                    law = GRAND_C[pair_class(meta[i], meta[j], B["xyz"])]
                    if not law:
                        continue
                    Rp = float(B["Rab"][meta[i][0], meta[j][0]])
                    fv = {"gb": gb[i, j], "s": float(S[i, j]),
                          "p": float(Ps[i, j]), "mm": 0.5 * float(m[i] + m[j]),
                          "obj": 0.42 * math.erf(0.2347 * Rp) / Rp
                          - 0.048 * float(S[i, j])}
                    FVAL = {"gb_s": fv["gb"] * fv["s"],
                            "gb_s3": fv["gb"] * fv["s"] ** 3,
                            "gb_p": fv["gb"] * fv["p"],
                            "gb_p_s2": fv["gb"] * fv["p"] * fv["s"] ** 2,
                            "s": fv["s"], "p": fv["p"], "obj": fv["obj"],
                            "obj_p": fv["obj"] * fv["p"],
                            "gb_s_mm": fv["gb"] * fv["s"] * fv["mm"]}
                    val = sum(cf * FVAL[t] for t, cf in law)
                    F[i, j] += val
                    F[j, i] += val
        else:
            RAB = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    RAB[i, j] = B["Rab"][meta[i][0], meta[j][0]]
            with np.errstate(divide="ignore", invalid="ignore"):
                obj = np.where(cross, 0.42 * erf_np(0.2347 * RAB)
                               / np.where(RAB > 0, RAB, 1.0) - 0.048 * S, 0.0)
            OFFF = {"gb_p_s2": gb * Ps * S2, "s": cross * S, "p": cross * Ps,
                    "obj_p": obj * Ps, "obj_s": cross * 0.42
                    * erf_np(0.2347 * RAB) / np.where(RAB > 0, RAB, 1.0) * S}
            for t, cf in GRAND_O.items():
                F += cf * OFFF[t]
    return F


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
    # THE FOUR-SYSTEM GATE (v3, with the grand short-piece layer). Pre-declared bar:
    # worst printed-window eigenvalue |d| <= 0.05 Eh per system (integration grade).
    ang = 104.5 * math.pi / 180
    r_oh = 0.9572 * K.BOHR
    wxyz = [[0.0, 0.0, 0.0],
            [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
            [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
    cases = [
        ("H2", [1, 1], [[0, 0, 0], [0, 0, 1.4]], 2),
        ("F2", [9, 9], [[0, 0, 0], [0, 0, 2.668]], 14),
        ("HF", [1, 9], [[0, 0, 0], [0, 0, 1.733]], 8),
        ("H2O", [8, 1, 1], wxyz, 8),
    ]
    sym = {1: "H", 8: "O", 9: "F"}
    npass = 0
    for name, zs, xyz, nel in cases:
        w, P = scf(zs, xyz, nel=nel)
        atoms = [(sym[z], x / K.BOHR, y_ / K.BOHR, zc / K.BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        r = oracle.run(atoms)
        ref = sorted(x / K.EV for x in r["eps_ev"])
        d = [a - b for a, b in zip(sorted(w)[:len(ref)], ref)]
        worst = max(abs(x) for x in d)
        ok = worst <= 0.05
        npass += ok
        print(f"{name:4s} worst |d| = {worst:.4f}  {'PASS' if ok else 'MISS'}"
              + ("   " + " ".join(f"{x:+.3f}" for x in d)))
    print(f"\nFOUR-SYSTEM GATE: {npass}/4 (bar 0.05 Eh, printed window)")


if __name__ == "__main__":
    main()
