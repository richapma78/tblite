"""f2_stretch.py -- the F2 stretch: multi-shell EHT decomposition on the reconstructed Fock.

Why F2 first (not HF): homonuclear means q = 0 exactly, so the ES off-diagonal shift that a
polar molecule adds to F is ZERO by symmetry -- the same trick that made H2 clean. And fluorine
brings two valence shells whose L4 differ hugely (s 0.7246 -> k = L4^2/2 = 0.2625; p 1.1269 ->
0.6350), so ONE molecule carries the per-AO vs per-pair question, the sigma/pi prefactor
question, and a second CN radius for the short-range candidates. HF comes after; once the F2
laws are fixed, HF's residual IS the ES off-diagonal measurement rather than a nuisance.

Decomposition per element (all from measurements gated elsewhere):
    F      <- fock_recon.fock_ao (gate: ortho 2.4e-9, exact pi zeros, banked-H2 tie 2.3e-9)
    X      <- Euler sum over the exchange slots: X_munu = sum_slots L5_slot * dF_munu/dL5_slot
              (degree-1 homogeneity, the H2-validated form, now two slots: L5(F,s), L5(F,p))
    A      <- the FULL analytic ACP matrix: A = V diag(c) V^T with V_i,g = <AO_i|g> computed by
              appending the L9 projector shells to the AO list and reading our gated overlap's
              cross block -- all onsite l-mismatch zeros arise by construction, no case-work
    EHT    =  F - X - A, charted per element type: (s,s') (s,pz') (pz,pz') (px,px') + the
              ONSITE (s,pz) null.

Pre-declared verdicts (written before the run):
  V-scale  On the CROSS element (Fs,Fpz'): 2-D scan of (k_s, k_p) in the scaled overlap
           S~(k_s,k_p); per-AO CONFIRMED iff the argmin lands within +-0.05 of
           (0.2625, 0.6350) AND off the k_s = k_p diagonal by > 2x the scan step; any per-pair
           rule would sit ON the diagonal. Same-shell elements can't discriminate (rules
           coincide); they get 1-D scans vs their own L4^2/2 (H2-style consistency).
  V-factor Geometric level-composition CONFIRMED iff T(s,p)^2/(T(ss)*T(pp)) = 1.00 +- 0.05
           on the tail (T = the flat EHT/S~ tail constant); arithmetic-flavored iff
           2T(s,p)/(T(ss)+T(pp)) hits 1.00 +- 0.05 with the geometric form failing; else
           report the two numbers and claim nothing.
  V-null   EHT(Fs,Fpz) ONSITE must stay < 5% of the two-center (Fs,Fpz') at r_e; a violation
           is recorded as a real onsite-mixing channel, not "fixed".
  V-pi     K_pi/K_sigma = T(px,px')/T(pz,pz') is REPORTED (no prior; Hueckel folklore says
           pi < sigma but folklore is not a gate).
Points where the SCF fails or the printed gap < 0.5 eV are excluded and listed (F2 at long R
is multireference; we do not use sick points).

Writes data/f2-stretch.json. Runs under the gpudft python in WSL.
"""
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
OUT = os.path.join(HERE, "data", "f2-stretch.json")

RS = [2.0, 2.3, 2.668, 3.0, 3.4, 3.8, 4.2, 4.6, 5.2, 6.0]      # Bohr
D_L5, D_C = 0.05, 0.02
# AO indices (labels gated: atom0 s,p(m0,m1,m2=z) = 0..3; atom1 = 4..7)
SS, SPZ, PZS, PZPZ, PXPX, ONSITE = (0, 4), (0, 7), (3, 4), (3, 7), (1, 5), (0, 3)


def f2(R_bohr):
    return [("F", 0, 0, 0), ("F", 0, 0, R_bohr / K.BOHR)]


def acp_matrix(zs, xyz, charge=0, parts=False):
    """A_munu = sum_projectors c * <mu|g><g|nu>, via projector shells appended to the AO list.
    parts=True additionally returns (V, per-column c, per-column (atom, l) labels)."""
    P = params.parse()
    shells, _ = overlap.build_shells(zs, xyz, charge=charge)
    n = sum(2 * s["l"] + 1 for s in shells)
    proj, cs, plab = [], [], []
    for at, z in enumerate(zs):
        l9 = P["element"][z]["l9"]
        for lp in range(4):
            c, zeta = l9[lp], l9[4 + lp]
            if abs(c) < 1e-12 or zeta <= 0:
                continue
            proj.append({"l": lp, "at": at, "exp": np.array([zeta]),
                         "coef": np.array([1.0])})
            cs += [c] * (2 * lp + 1)
            plab += [(at, lp)] * (2 * lp + 1)
    Sbig = overlap.overlap(list(zs), xyz, charge=charge, shells=shells + proj,
                           ao_order="pyscf")
    V = Sbig[:n, n:]
    A = V @ np.diag(cs) @ V.T
    return (A, V, cs, plab) if parts else A


def sweep():
    e = params.parse()["element"][9]
    L5s, L5p = e["shells"][3][0], e["shells"][3][1]
    rows, sick = [], []
    for R in RS:
        params.write_perturbed({})
        r0 = oracle.run(f2(R))
        if r0["gap_ev"] is None or r0["gap_ev"] < 0.5 or r0["scf_iterations"] is None:
            sick.append({"R": R, "gap_ev": r0["gap_ev"]})
            print(f"  R={R:4.2f}  EXCLUDED (gap {r0['gap_ev']} eV)", flush=True)
            continue
        base = fock_recon.fock_ao(f2(R))
        Fm = {}
        for tag, col, val in (("s+", 0, L5s + D_L5), ("s-", 0, L5s - D_L5),
                              ("p+", 1, L5p + D_L5), ("p-", 1, L5p - D_L5)):
            params.write_perturbed({(9, 5, col): val})
            Fm[tag] = fock_recon.fock_ao(f2(R))["F"]
        X = (L5s * (Fm["s+"] - Fm["s-"]) + L5p * (Fm["p+"] - Fm["p-"])) / (2 * D_L5)
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        A = acp_matrix([9, 9], xyz)
        EHT = base["F"] - X - A
        row = {"R": R, "gap_ev": r0["gap_ev"], "cn": r0["eeq"][0]["cn"],
               "S_ss": float(base["S"][SS]), "S_spz": float(base["S"][SPZ]),
               "S_pzpz": float(base["S"][PZPZ]), "S_pxpx": float(base["S"][PXPX]),
               "F_el": {k: float(base["F"][ix]) for k, ix in
                        (("ss", SS), ("spz", SPZ), ("pzs", PZS), ("pzpz", PZPZ),
                         ("pxpx", PXPX), ("onsite_spz", ONSITE))},
               "X_el": {k: float(X[ix]) for k, ix in
                        (("ss", SS), ("spz", SPZ), ("pzpz", PZPZ), ("pxpx", PXPX),
                         ("onsite_spz", ONSITE))},
               "A_el": {k: float(A[ix]) for k, ix in
                        (("ss", SS), ("spz", SPZ), ("pzpz", PZPZ), ("pxpx", PXPX),
                         ("onsite_spz", ONSITE))},
               "EHT_el": {k: float(EHT[ix]) for k, ix in
                          (("ss", SS), ("spz", SPZ), ("pzs", PZS), ("pzpz", PZPZ),
                           ("pxpx", PXPX), ("onsite_spz", ONSITE))}}
        rows.append(row)
        E = row["EHT_el"]
        print(f"  R={R:4.2f} gap {r0['gap_ev']:5.2f}  CN {row['cn']:.4f}  "
              f"EHT ss {E['ss']:+.5f}  spz {E['spz']:+.5f}  pzpz {E['pzpz']:+.5f}  "
              f"pxpx {E['pxpx']:+.5f}  onsite {E['onsite_spz']:+.5f}", flush=True)
    return rows, sick


def validate_acp(rs=(3.8, 2.668)):
    """c-slot FDs vs the analytic per-slot response dA/dc_l = sum_{cols of that l, BOTH atoms}
    V_col V_col^T (the slot is per-element). Tail must match; short range carries the known
    SCF-relaxation excess (measured on H2)."""
    e = params.parse()["element"][9]
    out, worst_tail = [], 0.0
    for R in rs:
        xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
        _A, V, _cs, plab = acp_matrix([9, 9], xyz, parts=True)
        rec = {"R": R}
        for lp, col in (("s", 0), ("p", 1)):
            c0 = e["l9"][col]
            params.write_perturbed({(9, 9, col): c0 + D_C})
            Fp = fock_recon.fock_ao(f2(R))["F"]
            params.write_perturbed({(9, 9, col): c0 - D_C})
            Fm = fock_recon.fock_ao(f2(R))["F"]
            dF = (Fp - Fm) / (2 * D_C)
            cols = [j for j, (_at, l) in enumerate(plab) if l == (0 if lp == "s" else 1)]
            dA = V[:, cols] @ V[:, cols].T
            rec[lp] = {}
            for k, ix in (("ss", SS), ("spz", SPZ), ("pzpz", PZPZ), ("pxpx", PXPX)):
                fd, an = float(dF[ix]), float(dA[ix])
                rec[lp][k] = {"fd": fd, "analytic": an}
                if R >= 3.8:
                    worst_tail = max(worst_tail, abs(fd - an))
            print(f"  R={R:4.2f} c_{lp}: " + "  ".join(
                f"{k} {rec[lp][k]['fd']:+.4f}/{rec[lp][k]['analytic']:+.4f}"
                for k in ("ss", "spz", "pzpz", "pxpx")), flush=True)
        out.append(rec)
    print(f"  ACP tail max|FD-analytic| = {worst_tail:.2e} "
          f"({'OK' if worst_tail < 5e-3 else 'MISMATCH -- decomposition suspect'})")
    return out, worst_tail


def scaled_S(R, kmap, charge=0):
    """Overlap with per-shell exponent scales: kmap[(at, l)] -> k."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    shells, _ = overlap.build_shells([9, 9], xyz, charge=charge)
    sc = [{"l": s["l"], "at": s["at"], "exp": s["exp"] * kmap[(s["at"], s["l"])],
           "coef": s["coef"].copy()} for s in shells]
    return overlap.overlap([9, 9], xyz, charge=charge, shells=sc, ao_order="pyscf")


def tail_T(rows, el, ix, ks, kp, tail):
    vals = []
    for r in rows:
        if r["R"] < tail:
            continue
        St = scaled_S(r["R"], {(0, 0): ks, (1, 0): ks, (0, 1): kp, (1, 1): kp})
        vals.append(r["EHT_el"][el] / float(St[ix]))
    m = float(np.mean(vals))
    spread = float((max(vals) - min(vals)) / abs(m)) if m else None
    return m, spread


def main():
    e = params.parse()["element"][9]
    k_s, k_p = e["shells"][2][0] ** 2 / 2, e["shells"][2][1] ** 2 / 2
    try:
        rows, sick = sweep()
        print("\nACP validation (c-slot FD vs analytic V V^T):", flush=True)
        acp_rec, acp_worst = validate_acp()
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored", flush=True)
    tail = 3.8
    live = [r for r in rows if r["R"] >= tail]
    print(f"\ntail = R >= {tail} ({len(live)} pts, CN there: "
          f"{['%.4f' % r['cn'] for r in live]})")

    # ---- internal consistency: (s,pz') vs (pz,s') antisymmetry under center swap
    anti = max(abs(r["EHT_el"]["spz"] + r["EHT_el"]["pzs"]) for r in rows)
    print(f"center-swap check: max|EHT(s,pz') + EHT(pz,s')| = {anti:.2e}")

    # ---- V-scale: 2-D scan on the cross element
    print("\nV-scale: 2-D (k_s, k_p) scan on EHT(s,pz'), tail flatness:")
    grid = np.arange(0.10, 0.96, 0.05)
    best = None
    for a in grid:
        for b in grid:
            vals = [r["EHT_el"]["spz"] /
                    float(scaled_S(r["R"], {(0, 0): a, (1, 0): a, (0, 1): b, (1, 1): b})[SPZ])
                    for r in live]
            m = np.mean(vals)
            sp = (max(vals) - min(vals)) / abs(m)
            if best is None or sp < best[2]:
                best = (a, b, sp, float(m))
    a0, b0, _, _ = best
    fine = np.arange(-0.04, 0.0401, 0.01)
    for a in a0 + fine:
        for b in b0 + fine:
            vals = [r["EHT_el"]["spz"] /
                    float(scaled_S(r["R"], {(0, 0): a, (1, 0): a, (0, 1): b, (1, 1): b})[SPZ])
                    for r in live]
            m = np.mean(vals)
            sp = (max(vals) - min(vals)) / abs(m)
            if sp < best[2]:
                best = (float(a), float(b), float(sp), float(m))
    ks_star, kp_star, sp_star, T_sp = best
    on_diag = abs(ks_star - kp_star) <= 0.02
    hit = abs(ks_star - k_s) <= 0.05 and abs(kp_star - k_p) <= 0.05
    verdict_scale = ("PER-AO CONFIRMED" if hit and not on_diag else
                     "ON-DIAGONAL (a per-pair rule)" if on_diag else
                     "NEITHER -- report raw")
    print(f"  argmin (k_s*, k_p*) = ({ks_star:.2f}, {kp_star:.2f})  spread {sp_star * 100:.1f}%  "
          f"T {T_sp:+.4f}\n  per-AO prediction ({k_s:.4f}, {k_p:.4f})  -> {verdict_scale}")

    # ---- same-shell 1-D consistency + T constants at the per-AO scales
    Ts = {}
    for el, ix, kk in (("ss", SS, (k_s, k_p)), ("pzpz", PZPZ, (k_s, k_p)),
                       ("pxpx", PXPX, (k_s, k_p)), ("spz", SPZ, (k_s, k_p))):
        T, sp = tail_T(rows, el, ix, kk[0], kk[1], tail)
        Ts[el] = T
        print(f"  T({el:5s}) at per-AO scales = {T:+.5f}  (tail spread {sp * 100:.1f}%)")

    # ---- V-factor
    geo = Ts["spz"] ** 2 / (Ts["ss"] * Ts["pzpz"])
    ari = 2 * Ts["spz"] / (Ts["ss"] + Ts["pzpz"])
    v_factor = ("GEOMETRIC" if abs(geo - 1) <= 0.05 else
                "ARITHMETIC-flavored" if abs(ari - 1) <= 0.05 else "NEITHER")
    print(f"\nV-factor: T(s,p)^2/(T(ss)T(pp)) = {geo:.3f}   2T(s,p)/(T(ss)+T(pp)) = {ari:.3f}"
          f"   -> {v_factor}")

    # ---- V-pi
    print(f"V-pi: K_pi/K_sigma = T(pxpx)/T(pzpz) = {Ts['pxpx'] / Ts['pzpz']:.3f}")

    # ---- V-null
    re_row = min(rows, key=lambda r: abs(r["R"] - 2.668))
    null_frac = abs(re_row["EHT_el"]["onsite_spz"]) / abs(re_row["EHT_el"]["spz"])
    print(f"V-null: onsite EHT(Fs,Fpz) at r_e = {re_row['EHT_el']['onsite_spz']:+.5f} = "
          f"{null_frac * 100:.1f}% of the two-center (bar: < 5%) "
          f"{'PASS' if null_frac < 0.05 else 'VIOLATED -- onsite mixing channel recorded'}")

    json.dump({"rows": rows, "sick": sick, "tail": tail,
               "acp_validation": {"records": acp_rec, "tail_worst": acp_worst},
               "scale_scan": {"ks_star": ks_star, "kp_star": kp_star, "spread": sp_star,
                              "per_ao_prediction": [k_s, k_p], "verdict": verdict_scale},
               "T": Ts, "factor": {"geometric": geo, "arithmetic": ari, "verdict": v_factor},
               "pi_over_sigma": Ts["pxpx"] / Ts["pzpz"],
               "onsite_null_frac": null_frac},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
