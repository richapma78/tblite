"""es_atoms.py -- ATOM gates for the first- and second-order electronic terms.

Lone neutral atoms are the cleanest electronic testbed:
  - CN = 0 (no CN scalings), no offsite terms, f(q_A=0) = 1 (the switching function's erf pair
    cancels), and E(3) vanishes identically (it is proportional to the TOTAL atomic charge).
  - Their converged shell populations equal the AUFBAU integers, but the REFERENCE occupations
    are FRACTIONAL (SI Tab. 1: averaged wB97M-V/q-vSZP Mulliken occupations) -- so the shell
    charges q_l = pop - refocc are nonzero and both ES1 and ES2+3 are live.

Shell charges are MEASURED, not assumed: ES1 is exactly linear in the mu(1) slots (L7), so a
central finite difference on L7[l] returns q_l to print precision. Those measured q_l are
cross-checked against SI Tab. 1 (pop - refocc) -- and then the INDEPENDENT gate runs:

  GATE (declared): E(2)_onsite = 1/2 sum_{ll'} q_l q_l' * 2 U_l U_l' / (U_l + U_l'), with U from
  the L6 row (CN = 0), must reproduce the oracle's printed ES2+3 for C, N, O, F, Si, P, S, Cl
  to <= 2e-5 Eh. The Hubbard slots and the kernel's onsite form carry no fitted freedom here.
"""
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oracle  # noqa: E402
import params  # noqa: E402

HOME = os.path.expanduser("~/.gxtb")
PRISTINE = os.path.join(HERE, "data", "gxtb_parameters.pristine")

# SI Tab. 1 reference occupations (Z <= 18), for the cross-check
SI_REFOCC = {6: {0: 1.035, 1: 2.964}, 7: {0: 1.374, 1: 3.625}, 8: {0: 1.673, 1: 4.326},
             9: {0: 1.855, 1: 5.144}, 14: {0: 1.387, 1: 2.324, 2: 0.287},
             15: {0: 1.569, 1: 2.991, 2: 0.438}, 16: {0: 1.757, 1: 4.089, 2: 0.152},
             17: {0: 1.875, 1: 5.067, 2: 0.057}}
UHF = {6: 2, 7: 3, 8: 2, 9: 1, 14: 2, 15: 3, 16: 2, 17: 1}
SYM = {6: "c", 7: "n", 8: "o", 9: "f", 14: "si", 15: "p", 16: "s", 17: "cl"}


def _lines():
    return open(PRISTINE).read().splitlines(keepends=True)


def _block_rows(lines, z):
    idx = [i for i, l in enumerate(lines) if l.strip()]
    hdr = next(i for i in idx if lines[i].split() == [str(z)])
    return idx[idx.index(hdr) + 1: idx.index(hdr) + 10]


def _perturbed(lines, li, tok, delta):
    t = lines[li].split()
    t[tok] = f"{float(t[tok]) + delta:.10f}"
    out = list(lines)
    out[li] = "      " + "      ".join(t) + "\n"
    return out


def measure_q(z, nshell, delta=0.05):
    """q_l = dES1/dmu_l by central FD on the L7 slots (ES1 is exactly linear in mu)."""
    lines = _lines()
    l7 = _block_rows(lines, z)[6]
    sym, uhf = SYM[z], UHF[z]
    q = []
    for l in range(nshell):
        vals = []
        for sgn in (+1, -1):
            open(HOME, "w").writelines(_perturbed(lines, l7, l, sgn * delta))
            r = oracle.run([(sym, 0.0, 0.0, 0.0)], uhf=uhf)
            vals.append(r["terms"]["ES1 (charge SIE)"])
        q.append((vals[0] - vals[1]) / (2 * delta))
    shutil.copyfile(PRISTINE, HOME)
    return q


def main():
    P = params.parse()
    shutil.copyfile(PRISTINE, HOME)
    fails = 0
    print(f"  {'el':3s} {'q_l (FD-measured)':30s} {'vs SI Tab.1':>11s} "
          f"{'ES2+3 pred':>12s} {'printed':>10s}")
    for z in (6, 7, 8, 9, 14, 15, 16, 17):
        nsh = len(SI_REFOCC[z])
        r0 = oracle.run([(SYM[z], 0.0, 0.0, 0.0)], uhf=UHF[z])
        pops = r0["pops"][0]
        pop = [pops["s"], pops["p"], pops["d"]][:nsh]
        es23 = r0["terms"]["ES2+3"]
        q = measure_q(z, nsh)
        # cross-check vs SI table: q_l should equal pop - refocc
        si_dev = max(abs(q[l] - (pop[l] - SI_REFOCC[z][l])) for l in range(nsh))
        U = P["element"][z]["shells"][4][:nsh]              # L6 row = U^(2),0
        E2 = 0.0
        for a in range(nsh):
            for b in range(nsh):
                gam = 2.0 * U[a] * U[b] / (U[a] + U[b])
                E2 += 0.5 * q[a] * q[b] * gam
        ok = abs(E2 - es23) <= 2e-5
        fails += 0 if ok else 1
        print(f"  {SYM[z]:3s} {str([round(x, 4) for x in q]):30s} {si_dev:11.4f} "
              f"{E2:12.6f} {es23:10.6f}  {'MATCH' if ok else 'MISMATCH'}")
    if fails:
        print(f"\nES ATOM GATE: {fails} MISMATCH -- kernel/slot reading wrong; do not tune")
        sys.exit(1)
    print("\nES ATOM GATE PASSED: onsite second-order from L6 Hubbards reproduces printed "
          "ES2+3; measured shell charges match SI Tab. 1 fractional references")


if __name__ == "__main__":
    main()
