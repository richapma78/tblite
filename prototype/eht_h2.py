"""eht_h2.py -- the H2 off-diagonal decomposition, COMMITTED, plus the three probes.

Push 34 extracted the naked EHT core from H2 in a one-off that never landed in the repo (only
its unlabeled table, data/h2-h0.json, did). This script IS that instrument, rebuilt with labeled
output, and it carries the three follow-up probes:

  PROBE 1 (fine-scan): the Hamiltonian-basis exponent scale. Build S~12(k) = overlap with all
    primitive exponents x k (adapted coefficients kept, everything renormalized) and find the k
    that makes EHT12/S~12(k) FLATTEST over the tail (R >= 2.5 Bohr). Push 34 saw ~0.6 at +-1%.
  PROBE 2 (FD the hypothesis): 0.6 ~ L4(H)^2/2 = 0.59846 was flagged HYPOTHESIS ONLY (the G1[2]
    numerology lesson). Perturb L4(H) by +0.10 and L1[7](H) by x1.1 IN THE PARAMETER FILE,
    re-extract the core, re-scan k*. Pre-declared verdicts, written before the run:
      - L4 hypothesis CONFIRMED  iff baseline k* in [0.586, 0.611] AND k*(L4+0.1) in
        [0.693, 0.733] (prediction (1.1940^2)/2 = 0.71283).
      - L4 hypothesis REFUTED    iff k*(L4+0.1) within +-0.02 of baseline while the plateau
        amplitude moves by > 3%.
      - L1[7]: amplitude-only iff k* fixed and the plateau scales; scale iff k* moves.
        NOTE (post-run correction, kept for the record): this arm's first draft carried an
        "L1[7]^2-scale" band [0.782, 0.822] built on a MISREAD slot value (0.8142 -- that is
        L1[6] = kU; the true L1[7] is 0.0224), so that band was void before the run. The
        discriminator actually used -- k* moved vs pinned -- is value-independent and stands.
      - anything else: report the raw movement, claim nothing.
  PROBE 3 (short-range ingredient): chart ratio(R) = EHT12/S~12(k*) against the CN channels the
    method already carries (the printed Eq.47 CN and the basis CN come free with every oracle
    run). A candidate WINS only if its 1-parameter fit rms beats the runner-up by 2x and the
    residual is at the tail's own +-1% noise; otherwise report the table and claim nothing.

The decomposition (per R, all measured):
    F11, F12   <- printed eigenvalue pair + our machine-exact S12 (2x2 inversion; H is one
                  contracted s AO per atom in q-vSZP, checked)
    X12        <- L5-linearity: X12 = L5 * dF12/dL5 (central FD on the file slot)
    ACP12      <- ANALYTIC two-center: sum over centers C of c_s*<AO1|g_C><g_C|AO2> with
                  g = normalized s Gaussian(zeta_s from L9); VALIDATED per R against the c_s FD
                  response before it is subtracted (push 34's single-point check, curve-wide).
                  The analytic value is the BARE matrix element -- at R < 2 the FD response
                  exceeds it by up to ~4% (SCF density relaxation, measured, recorded); the
                  decomposition subtracts the bare element by construction.
    EHT12      = F12 - X12 - ACP12   -- the naked Hueckel core

Runs in WSL under the gpudft python (scipy for the c2s builder). Writes
data/h2-eht-probes.json (labeled). Modes: --check (R=2.5 selftest), default (baseline sweep +
probes 1+3), --probe2 (the two perturbation sweeps; reads the baseline JSON).
"""
import copy
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as K  # noqa: E402
import oracle  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
OUT = os.path.join(HERE, "data", "h2-eht-probes.json")

RS = [0.9, 1.0, 1.2, 1.4, 1.7, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]   # Bohr
TAIL = [r for r in RS if r >= 2.5]
RS_P2 = [1.4, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
D_L5, D_CS = 0.05, 0.02


# --------------------------------------------------------------- parameter-file slot editor
def _write_edits(edits):
    """edits: {(row_offset_from_header, col): value} on hydrogen's block; {} -> pristine."""
    lines = open(PRISTINE).read().splitlines(keepends=True)
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == ["1"])
    for (off, col), val in edits.items():
        li = idx[idx.index(hdr) + off]
        t = lines[li].split()
        t[col] = f"{val:.10f}"
        lines[li] = "      " + "      ".join(t) + "\n"
    open(HOME, "w").writelines(lines)


# --------------------------------------------------------------------------- the inversion
def fock(R, edits=None):
    """One oracle run at separation R (Bohr) -> (e_bond, e_anti, S12, F11, F12, cn, cn_basis)."""
    if edits is not None:
        _write_edits(edits)
    r = oracle.run([("H", 0, 0, 0), ("H", 0, 0, R / K.BOHR)])
    eps = [e / K.EV for e in r["eps_ev"]]
    assert len(eps) == 2, f"expected 2 eigenvalues, got {len(eps)}"
    eb, ea = min(eps), max(eps)
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    S = overlap.overlap([1, 1], xyz)[0, 1]
    F11 = ((1 + S) * eb + (1 - S) * ea) / 2
    F12 = ((1 + S) * eb - (1 - S) * ea) / 2
    return {"e_bond": eb, "e_anti": ea, "S12": float(S), "F11": F11, "F12": F12,
            "cn": r["eeq"][0]["cn"], "cn_basis": r["eeq"][0]["cn_basis"]}


# ------------------------------------------------------------------- analytic ACP integrals
def acp_integrals(R):
    """(onsite, cross) projector overlaps <AO1|g_A>, <AO1|g_B> for the ADAPTED H2 basis at R;
    g = normalized s Gaussian with zeta_s from L9. Both-center ACP responses follow:
        dF12/dc_s = 2*on*tw        dF11/dc_s = on^2 + tw^2"""
    e = params.parse()["element"][1]
    zeta = e["l9"][4]
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    shells, _ = overlap.build_shells([1, 1], xyz)
    sh = shells[0]                                       # AO1: the s shell on atom A
    c = sh["coef"] * np.array([overlap.prim_norm(a, 0) for a in sh["exp"]])
    s_rad = (np.pi ** 1.5) / (sh["exp"][:, None] + sh["exp"][None, :]) ** 1.5
    c = c / math.sqrt(float(c @ s_rad @ c))
    ng = overlap.prim_norm(zeta, 0)
    on = sum(ca * ng * overlap.cart_prims(a, zeta, xyz[0], xyz[0], 0, 0)[0, 0]
             for ca, a in zip(c, sh["exp"]))
    tw = sum(ca * ng * overlap.cart_prims(a, zeta, xyz[0], xyz[1], 0, 0)[0, 0]
             for ca, a in zip(c, sh["exp"]))
    return float(on), float(tw)


# ------------------------------------------------------------------------ scaled overlaps
def s_scaled(R, k):
    """S~12(k): every primitive exponent x k, adapted coefficients kept, renormalized."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    shells, _ = overlap.build_shells([1, 1], xyz)
    scaled = [{"l": s["l"], "at": s["at"], "exp": s["exp"] * k, "coef": s["coef"].copy()}
              for s in shells]
    return float(overlap.overlap([1, 1], xyz, shells=scaled)[0, 1])


def scan_scale(rows, ks):
    """Flatness of EHT12/S~12(k) over the tail -> [(k, spread, mean_ratio)]."""
    out = []
    tail = [r for r in rows if r["R"] >= 2.5]
    for k in ks:
        ratios = [r["EHT12"] / s_scaled(r["R"], k) for r in tail]
        m = float(np.mean(ratios))
        out.append((float(k), float((max(ratios) - min(ratios)) / abs(m)), m))
    return out


def best_scale(rows):
    coarse = scan_scale(rows, np.arange(0.40, 1.001, 0.02))
    k0 = min(coarse, key=lambda t: t[1])[0]
    fine = scan_scale(rows, np.arange(k0 - 0.03, k0 + 0.0301, 0.002))
    return min(fine, key=lambda t: t[1]), coarse, fine


# ------------------------------------------------------------------------------- the sweep
def sweep(rs, base_edits=None, validate_acp=True, label="baseline"):
    e = params.parse()["element"][1]
    L5, CS = e["shells"][3][0], e["l9"][0]
    if base_edits:                                       # probe-2: the slot itself moves
        L5 = base_edits.get((5, 0), L5)
    rows = []
    for R in rs:
        base = fock(R, edits=dict(base_edits or {}))
        p = fock(R, edits={**(base_edits or {}), (5, 0): L5 + D_L5})
        m = fock(R, edits={**(base_edits or {}), (5, 0): L5 - D_L5})
        dF12_dL5 = (p["F12"] - m["F12"]) / (2 * D_L5)
        X12 = L5 * dF12_dL5
        on, tw = acp_integrals(R)
        resp12, resp11 = 2 * on * tw, on * on + tw * tw
        row = {"R": R, **base, "dF12_dL5": dF12_dL5, "X12": X12,
               "acp_on": on, "acp_tw": tw, "acp_resp12_analytic": resp12}
        if validate_acp:
            pc = fock(R, edits={**(base_edits or {}), (9, 0): CS + D_CS})
            mc = fock(R, edits={**(base_edits or {}), (9, 0): CS - D_CS})
            row["dF12_dcs"] = (pc["F12"] - mc["F12"]) / (2 * D_CS)
            row["dF11_dcs"] = (pc["F11"] - mc["F11"]) / (2 * D_CS)
            row["acp_resp11_analytic"] = resp11
        row["ACP12"] = CS * resp12
        row["EHT12"] = row["F12"] - X12 - row["ACP12"]
        rows.append(row)
        v = (f"  fd12 {row['dF12_dcs']:+.4f}/an {resp12:+.4f}  fd11 {row['dF11_dcs']:+.4f}"
             f"/an {resp11:+.4f}") if validate_acp else ""
        print(f"  [{label}] R={R:4.2f}  S {base['S12']:+.5f}  F12 {base['F12']:+.6f}  "
              f"X12 {X12:+.6f}  ACP12 {row['ACP12']:+.6f}  EHT12 {row['EHT12']:+.6f}{v}",
              flush=True)
    return rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    e = params.parse()["element"][1]
    L4, L17 = e["shells"][2][0], e["l1"][7]
    try:
        if mode == "--check":
            _write_edits({})
            r = fock(2.5)
            ref = [0.4072334609268092, -0.2803005270651893, -0.24940460789490132]
            print(f"R=2.5: S12 {r['S12']:.10f} (banked {ref[0]:.10f})  "
                  f"F11 {r['F11']:+.7f} ({ref[1]:+.7f})  F12 {r['F12']:+.7f} ({ref[2]:+.7f})")
            assert abs(r["S12"] - ref[0]) < 1e-9 and abs(r["F11"] - ref[1]) < 2e-5 \
                and abs(r["F12"] - ref[2]) < 2e-5
            row = sweep([2.5], validate_acp=True, label="check")[0]
            da = abs(row["dF12_dcs"] - row["acp_resp12_analytic"])
            db = abs(row["dF11_dcs"] - row["acp_resp11_analytic"])
            print(f"ACP analytic-vs-FD: off-diag {da:.2e}  diag {db:.2e}")
            # L5 linearity: half/double delta
            L5 = e["shells"][3][0]
            f1 = fock(2.5, edits={(5, 0): L5 + 2 * D_L5})
            f2 = fock(2.5, edits={(5, 0): L5 - 2 * D_L5})
            wide = (f1["F12"] - f2["F12"]) / (4 * D_L5)
            print(f"L5 slope delta={D_L5}: {row['dF12_dL5']:+.6f}  2delta: {wide:+.6f}  "
                  f"curvature {abs(wide - row['dF12_dL5']):.2e}")
            print("CHECK PASSED" if da < 5e-3 and db < 5e-3 else "CHECK: ACP mismatch")
            return

        if mode == "--probe2":
            data = json.load(open(OUT))
            base_k, base_plateau = data["probe1"]["k_star"], data["probe1"]["plateau"]
            res = {}
            for tag, edits, pred in (
                    ("L4+0.10", {(4, 0): L4 + 0.10}, (L4 + 0.10) ** 2 / 2),
                    ("L1[7]x1.1", {(1, 7): L17 * 1.1}, None)):
                rows = sweep(RS_P2, base_edits=edits, validate_acp=False, label=tag)
                (k, spread, plat), _, _ = best_scale(rows)
                ptxt = f"scale-hypothesis prediction {pred:.3f}" if pred else \
                    "discriminator: k* moved (scale) vs pinned (amplitude)"
                print(f"  [{tag}] k* = {k:.3f} (spread {spread:.3f})  plateau {plat:+.4f}   {ptxt}")
                res[tag] = {"rows": rows, "k_star": k, "spread": spread, "plateau": plat,
                            "prediction_if_scale": pred}
            # pre-declared verdicts
            kL4 = res["L4+0.10"]["k_star"]
            if 0.586 <= base_k <= 0.611 and 0.693 <= kL4 <= 0.733:
                v4 = "CONFIRMED: k* moved to (L4+d)^2/2"
            elif abs(kL4 - base_k) <= 0.02 and \
                    abs(res["L4+0.10"]["plateau"] / base_plateau - 1) > 0.03:
                v4 = "REFUTED: k* fixed, amplitude moved -- L4 is a prefactor"
            else:
                v4 = "NEITHER: raw movement reported, no claim"
            k17 = res["L1[7]x1.1"]["k_star"]
            if abs(k17 - base_k) <= 0.02:
                ratio = res["L1[7]x1.1"]["plateau"] / base_plateau
                v17 = f"amplitude-only: plateau x{ratio:.3f} for a +10% slot move (k* pinned)"
            else:
                v17 = f"k* moved to {k17:.3f}: a scale participant -- characterize before claiming"
            print(f"\nVERDICT L4:    {v4}\nVERDICT L1[7]: {v17}")
            data["probe2"] = {**res, "verdict_L4": v4, "verdict_L17": v17}
            json.dump(data, open(OUT, "w"), indent=1)
            return

        # ------------------------------------------------------------ baseline + probes 1, 3
        rows = sweep(RS, validate_acp=True, label="baseline")
        acp_tail = max(abs(r["dF12_dcs"] - r["acp_resp12_analytic"]) for r in rows if r["R"] >= 2.0)
        acp_short = max((r["dF12_dcs"] - r["acp_resp12_analytic"]) / r["acp_resp12_analytic"]
                        for r in rows if r["R"] < 2.0)
        print(f"\nACP form validation (tail, R>=2.0): max |FD - analytic| = {acp_tail:.2e} "
              f"({'OK' if acp_tail < 5e-3 else 'MISMATCH -- do not subtract'})")
        print(f"ACP short-range FD excess over bare matrix element: up to {acp_short * 100:+.1f}% "
              f"(the SCF density-relaxation channel; the decomposition subtracts the BARE element "
              f"by construction, so probe-3's short-range chart carries this as a known systematic)")

        (k, spread, plateau), coarse, fine = best_scale(rows)
        hyp = L4 ** 2 / 2
        print(f"\nPROBE 1: k* = {k:.3f}  spread {spread * 100:.2f}%  plateau {plateau:+.4f}")
        print(f"  L4^2/2 = {hyp:.5f} -> {'within' if abs(k - hyp) <= 0.01 else 'OUTSIDE'} "
              f"+-0.01 of k* (consistency only; probe 2 decides)")

        print("\nPROBE 3: ratio(R) = EHT12/S~12(k*) vs 1-parameter candidates "
              "ratio = p*(1 + a*x):")
        ratios = np.array([r["EHT12"] / s_scaled(r["R"], k) for r in rows])
        cands = {"CN": np.array([r["cn"] for r in rows]),
                 "sqrt(CN)": np.array([math.sqrt(max(r["cn"], 0)) for r in rows]),
                 "CN_basis": np.array([r["cn_basis"] for r in rows]),
                 "sqrt(CN_basis)": np.array([math.sqrt(max(r["cn_basis"], 0)) for r in rows]),
                 "S12": np.array([r["S12"] for r in rows]),
                 "1/R": np.array([1 / r["R"] for r in rows])}
        table = {}
        for name, x in cands.items():
            A = np.vstack([np.ones_like(x), x]).T
            (p, pa), *_ = np.linalg.lstsq(A, ratios, rcond=None)
            resid = A @ np.array([p, pa]) - ratios
            table[name] = {"p": float(p), "a": float(pa / p), "rms": float(np.sqrt(np.mean(resid ** 2))),
                           "max": float(np.abs(resid).max())}
            print(f"  {name:15s} p {p:+.4f}  a {pa / p:+.4f}  rms {table[name]['rms']:.5f}  "
                  f"max {table[name]['max']:.5f}")
        ranked = sorted(table.items(), key=lambda t: t[1]["rms"])
        win, second = ranked[0], ranked[1]
        tail_noise = spread * abs(plateau)
        decided = win[1]["rms"] * 2 <= second[1]["rms"] and win[1]["max"] <= 2 * tail_noise
        print(f"  -> {'WINNER: ' + win[0] if decided else 'NO WINNER (pre-declared bar not met)'}"
              f"   [bar: rms 2x under runner-up AND max resid <= {2 * tail_noise:.4f}]")

        json.dump({"rows": rows,
                   "probe1": {"k_star": k, "spread": spread, "plateau": plateau,
                              "hypothesis_L4sq_over_2": hyp, "coarse": coarse, "fine": fine},
                   "probe3": {"ratios": ratios.tolist(), "fits": table,
                              "winner": win[0] if decided else None}},
                  open(OUT, "w"), indent=1)
        print(f"\nwrote {OUT}")
    finally:
        shutil.copyfile(PRISTINE, HOME)
        print("param file restored")


if __name__ == "__main__":
    main()
