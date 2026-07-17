"""hf_es_probe.py -- measure WHICH charges drive the ES terms on a polar system.

ES1 = sum_A sum_l mu_l * f(q_A) * q_l is exactly mu-linear, so the slope of the PRINTED
'ES1 (charge SIE)' against each mu slot reads the operative f(q)*q_l per shell directly:

    dES1/dmu_l(A) = f(q_A) * q_l(A)

Candidates on HF at r_e (3.5x apart -- decisive): EEQ atom charges (+-0.116) under some
shell partition, vs Mulliken shell charges (q_H = +0.41; q_Fs, q_Fp). The same round FDs
'ES2+3' against the U slots (the kernel's q^2 structure) as the cross-check, and the
difference between the measured ES2+3 and our onsite model with the measured charges
isolates the OFFSITE gamma2 contribution (target 3) as a bonus. Writes data/hf-es-probe.json.
"""
import json
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
D = 0.02
SLOTS = [("mu_Hs", (1, 7, 0)), ("mu_Fs", (9, 7, 0)), ("mu_Fp", (9, 7, 1)),
         ("U_Hs", (1, 6, 0)), ("U_Fs", (9, 6, 0)), ("U_Fp", (9, 6, 1))]


def term(raw, name):
    m = re.search(rf"^\s*{re.escape(name)}\s*:\s*(-?\d+\.\d+)", raw, re.M)
    return float(m.group(1))


def main():
    P = params.parse()
    vals = {"mu_Hs": P["element"][1]["shells"][5][0],
            "mu_Fs": P["element"][9]["shells"][5][0],
            "mu_Fp": P["element"][9]["shells"][5][1],
            "U_Hs": P["element"][1]["shells"][4][0],
            "U_Fs": P["element"][9]["shells"][4][0],
            "U_Fp": P["element"][9]["shells"][4][1]}
    out = {}
    try:
        for R in (1.733, 2.7):
            atoms = [("H", 0, 0, 0), ("F", 0, 0, R / K.BOHR)]
            params.write_perturbed({})
            r0 = oracle.run(atoms)
            eeq = [e["q"] for e in r0["eeq"]]
            rec = {"eeq": eeq, "slopes": {}}
            for name, key in SLOTS:
                params.write_perturbed({key: vals[name] + D})
                rp = oracle.run(atoms)
                params.write_perturbed({key: vals[name] - D})
                rm = oracle.run(atoms)
                tname = "ES1 (charge SIE)" if name.startswith("mu") else "ES2+3"
                sl = (term(rp["raw"], tname) - term(rm["raw"], tname)) / (2 * D)
                rec["slopes"][name] = sl
                print(f"  R={R} d{'ES1' if name.startswith('mu') else 'ES23'}/d{name} = "
                      f"{sl:+.5f}", flush=True)
            out[str(R)] = rec
            print(f"  EEQ q: H {eeq[0]:+.4f}  F {eeq[1]:+.4f}")
            print(f"  => operative f*q per shell: Hs {rec['slopes']['mu_Hs']:+.4f}  "
                  f"Fs {rec['slopes']['mu_Fs']:+.4f}  Fp {rec['slopes']['mu_Fp']:+.4f}")
            print(f"     candidates: EEQ-atom H {eeq[0]:+.4f}; Mulliken-shell "
                  f"(from the energy-gate run: H +0.41; F split)", flush=True)
    finally:
        shutil.copyfile(PRISTINE, os.path.expanduser("~/.gxtb"))
        print("param file restored")
    json.dump(out, open(os.path.join(HERE, "data", "hf-es-probe.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
