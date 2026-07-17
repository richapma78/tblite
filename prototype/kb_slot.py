"""kb_slot.py -- locate the Ham-basis sqrt(CN) adaptation knob (kb = 0.17 measured on H2).

The candidates (unnamed slots with CN-family or H2-visible responses in the atlas):
G1[4], L8[3](H), L8[7](H), and L8[2](H) -- the last one to EXECUTE a numerology
(L8[2](H)*b_density(H) = 0.720*0.2272 = 0.164 ~ kb; the G1[2] rule says test, never trust).

Method: central FD of F12(H2) per slot at R = 0.9, 1.4, 2.0 (basis-CN alive), compared to the
ANALYTIC dF12/dkb = amp*(-L2)*Pi*(dS~/dkb) at the fitted metric (k = 1.110, kb = 0.17,
b = 0.0315, amp = kW_s*kdiat_s(H)). Verdict (pre-declared): slot X IS the knob iff
resp_X(R)/analytic(R) is R-constant within +-10% over the three points; the constant is then
dkb/dX. If the constant ~ b_density(H) = 0.2272, the slot SCALES the density-basis b (the
SI's 'we scale ... the three parameters'). If no candidate matches, kb joins the hardcoded
ledger (like the k^shp,l globals). Writes data/kb-slot.json.
"""
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import adapt  # noqa: E402
import basisq  # noqa: E402
import fock_recon  # noqa: E402
import overlap  # noqa: E402
import params  # noqa: E402

BOHR = 1.8897261254578281
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
RS = (0.9, 1.4, 2.0)
D = 0.02
K_FIT, KB_FIT, B_SHP = 1.110, 0.17, 0.0315

P = params.parse()
eH = P["element"][1]
AMP = P["globals"]["g1"][0] * eH["l8"][0]          # kW_s * kdiat_sigma(H)
L2 = eH["shells"][0][0]

Bq = basisq.parse()
prims = Bq[1]["shells"][0][1]
exps0 = np.array([p[0] for p in prims])
c0v = np.array([p[1] for p in prims])
c1v = np.array([p[2] for p in prims])


def s_ham(R, kb):
    xyz = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, R]])
    cn = adapt.basis_cn([1, 1], xyz)[0]
    coef = c0v + c1v * (kb * math.sqrt(cn))
    sh = [{"l": 0, "at": i, "exp": exps0 * K_FIT, "coef": coef.copy()} for i in (0, 1)]
    return float(overlap.overlap([1, 1], xyz, shells=sh)[0, 1])


def analytic_dkb(R):
    return AMP * (-L2) * (1 + B_SHP * R) ** 2 * (s_ham(R, KB_FIT + 0.01) -
                                                 s_ham(R, KB_FIT - 0.01)) / 0.02


CANDS = [("G1[4]", ("g", 1, 4), P["globals"]["g1"][4]),
         ("L8[2]H", (1, 8, 2), eH["l8"][2]),
         ("L8[3]H", (1, 8, 3), eH["l8"][3]),
         ("L8[7]H", (1, 8, 7), eH["l8"][7])]


def main():
    ana = {R: analytic_dkb(R) for R in RS}
    print("analytic dF12/dkb at the fitted metric: " +
          "  ".join(f"R={R}: {ana[R]:+.4f}" for R in RS))
    out = {"analytic": {str(R): ana[R] for R in RS}, "slots": {}}
    try:
        for name, key, val in CANDS:
            row = {}
            for R in RS:
                atoms = [("H", 0, 0, 0), ("H", 0, 0, R / BOHR)]
                params.write_perturbed({key: val + D})
                Fp = fock_recon.fock_ao(atoms)["F"][0, 1]
                params.write_perturbed({key: val - D})
                Fm = fock_recon.fock_ao(atoms)["F"][0, 1]
                row[R] = float((Fp - Fm) / (2 * D))
            consts = [row[R] / ana[R] for R in RS]
            spread = (max(consts) - min(consts)) / abs(np.mean(consts)) \
                if abs(np.mean(consts)) > 1e-4 else None
            verdict = ("KNOB (dkb/dslot = %.4f)" % np.mean(consts)
                       if spread is not None and spread <= 0.10 and
                       abs(np.mean(consts)) > 0.01 else "not the knob")
            print(f"  {name:7s} resp " + "  ".join(f"{row[R]:+8.5f}" for R in RS) +
                  "   /analytic " + "  ".join(f"{c:+7.3f}" for c in consts) +
                  f"   [{verdict}]", flush=True)
            out["slots"][name] = {"resp": {str(R): row[R] for R in RS},
                                  "ratio": consts, "verdict": verdict}
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored")
    json.dump(out, open(os.path.join(HERE, "data", "kb-slot.json"), "w"), indent=1)
    print("wrote data/kb-slot.json")


if __name__ == "__main__":
    main()
