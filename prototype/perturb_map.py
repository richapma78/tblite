"""perturb_map.py -- name the gxtb_parameters slots by PERTURBING them one at a time.

The oracle reads ~/.gxtb at every run and prints every energy term, the EEQ charges, the AO
adaptation, eigenvalues and the gap. So: copy the parameter file, nudge ONE float by a small
delta, run fixed probes, and record WHICH observables moved. That maps file slots onto physical
roles mechanically -- no fitting, no guessing. The original file is restored afterwards (and on
any crash: the write happens fresh from the pristine copy each iteration).

Probes: H2 @ 1.4 bohr (repulsion/EHT on the simplest pair), HCl @ 2.409 bohr (polar pair:
charge-dependent terms light up), lone O atom (pure atomic terms).

Observables tracked per run: all printed terms, first/last eigenvalue, gap, q(atom 1),
CN(atom 1), q_eff total (atom 1), SCF iteration count.
"""
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oracle  # noqa: E402

GXTB_HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")
BOHR = 1.8897261254578281

PROBES = {
    "H2": ([("H", 0, 0, 0), ("H", 0, 0, 1.4 / BOHR)], 0, 0),
    "H2far": ([("H", 0, 0, 0), ("H", 0, 0, 4.0 / BOHR)], 0, 0),      # CN ~ 0: kCN goes dark here
    "HCl": ([("Cl", 0, 0, 0), ("H", 0, 0, 2.409 / BOHR)], 0, 0),
    "NH4+": ([("N", 0, 0, 0), ("H", 1.0134 / 1.0, 0, 0), ("H", -0.3378, 0.9553, 0),
              ("H", -0.3378, -0.4776, -0.8273), ("H", -0.3378, -0.4776, 0.8273)], 1, 0),
    "O": ([("O", 0, 0, 0)], 0, 2),
}


def observe():
    out = {}
    for name, (atoms, chg, uhf) in PROBES.items():
        r = oracle.run(atoms, charge=chg, uhf=uhf)
        o = {("t:" + k): v for k, v in r["terms"].items()}
        o["total"] = r["total"]
        o["eps0"] = r["eps_ev"][0] if r["eps_ev"] else None
        o["gap"] = r["gap_ev"]
        if r["eeq"]:
            o["q1"] = r["eeq"][0]["q"]
            o["cn1"] = r["eeq"][0]["cn"]
            o["cnb1"] = r["eeq"][0]["cn_basis"]
        if r["ao_adapt"]:
            o["qeff1"] = r["ao_adapt"][0]["total"]
        out[name] = o
    return out


def load_lines():
    return open(PRISTINE).read().splitlines(keepends=True)


def write_perturbed(lines, line_idx, tok_idx, delta):
    toks = lines[line_idx].split()
    val = float(toks[tok_idx]) + delta
    toks[tok_idx] = f"{val:.10f}"
    out = list(lines)
    out[line_idx] = "      " + "      ".join(toks) + "\n"
    with open(GXTB_HOME, "w") as f:
        f.writelines(out)


def diff(base, new, tol=1e-9):
    hits = []
    for probe in base:
        for k, v0 in base[probe].items():
            v1 = new[probe].get(k)
            if v0 is None or v1 is None:
                continue
            if abs(v1 - v0) > tol * max(1.0, abs(v0)):
                hits.append(f"{probe}.{k}({v1 - v0:+.2e})")
    return hits


def main():
    os.makedirs(os.path.dirname(PRISTINE), exist_ok=True)
    if not os.path.exists(PRISTINE):
        shutil.copyfile(GXTB_HOME, PRISTINE)
    lines = load_lines()

    # locate structural anchors: globals = first two non-blank lines; H block follows the lone "1"
    idx_nonblank = [i for i, l in enumerate(lines) if l.strip()]
    g1, g2 = idx_nonblank[0], idx_nonblank[1]
    h_hdr = next(i for i in idx_nonblank if lines[i].split() == ["1"])
    h_rows = idx_nonblank[idx_nonblank.index(h_hdr) + 1: idx_nonblank.index(h_hdr) + 10]

    targets = []
    for j in range(10):
        targets.append((f"G1[{j}]", g1, j))
    for j in range(10):
        targets.append((f"G2[{j}]", g2, j))
    for j in range(10):
        targets.append((f"H.L1[{j}]", h_rows[0], j))
    for li in range(1, 7):
        targets.append((f"H.L{li+1}[0]", h_rows[li], 0))
    for j in range(8):
        targets.append((f"H.L8[{j}]", h_rows[7], j))
    for j in range(8):
        targets.append((f"H.L9[{j}]", h_rows[8], j))

    shutil.copyfile(PRISTINE, GXTB_HOME)
    base = observe()
    results = {}
    try:
        for name, li, tj in targets:
            write_perturbed(lines, li, tj, 0.02)
            try:
                new = observe()
                results[name] = diff(base, new)
            except Exception as e:
                results[name] = [f"RUN-FAILED: {e}"]
            print(f"  {name:10s} -> " + ("; ".join(results[name][:6]) if results[name]
                                          else "(no change)"), flush=True)
    finally:
        shutil.copyfile(PRISTINE, GXTB_HOME)
        print("\n  pristine parameter file RESTORED")
    json.dump(results, open(os.path.join(HERE, "data", "perturb_map.json"), "w"), indent=1)
    print("  wrote data/perturb_map.json")


if __name__ == "__main__":
    main()
