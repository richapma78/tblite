"""restart.py -- parse the oracle's gxtbrestart file: the CONVERGED state as an input.

Why this matters: with the converged density matrix P in hand, every electronic term (Mulliken
populations, ES1, ES2+3, multipole ES, the Mulliken-approximated exchange) can be gated
INDEPENDENTLY against the oracle's printed per-term decomposition -- before any SCF loop of our
own exists. The SCF then becomes fixed-point iteration around already-gated pieces.

Format (decoded from the water file, 472 bytes, and pinned by the gate below):
  Fortran unformatted sequential records (4-byte length markers):
  record 1: nsao*(nsao+1)/2 doubles -- the symmetric DENSITY matrix, packed triangular
  record 2: nsao*nsao       doubles -- a full square matrix (MO coefficients)
  (open-shell files may carry more records -- probed when needed)

GATE (declared): Mulliken shell populations q_l = sum_{mu in l} (P S)_mu,mu computed with the
parsed P and OUR overlap engine's S must reproduce the oracle's printed per-shell populations
('a+b dens' block) for water AND HCl AND AcCl -- which simultaneously pins the packing order,
the AO ordering, and proves the oracle's AO normalization convention equals ours.
"""
import os
import struct
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

GXTB_V1 = os.environ.get("GXTB_V1_EXE", "/opt/gxtb-v1/gxtb")
BOHR = 1.8897261254578281


def read_records(path):
    """All Fortran sequential records as raw byte strings."""
    out = []
    with open(path, "rb") as f:
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (n,) = struct.unpack("<i", head)
            payload = f.read(n)
            f.read(4)                                   # trailer
            out.append(payload)
    return out


def unpack_symmetric(vec, n, order="row_lower"):
    """Packed triangle -> full symmetric matrix. row_lower: (1,1),(2,1),(2,2),(3,1)..."""
    M = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            M[i, j] = M[j, i] = vec[k]
            k += 1
    return M


def converged_state(atoms, charge=0, uhf=0):
    """Run the oracle in a scratch dir, parse its restart: dict(P, C, nsao, raw records)."""
    import oracle
    with tempfile.TemporaryDirectory() as d:
        lines = ["$coord"]
        for sym, x, y, z in atoms:
            lines.append(f" {x * BOHR:.10f} {y * BOHR:.10f} {z * BOHR:.10f} {sym.lower()}")
        lines.append("$end")
        open(os.path.join(d, "coord"), "w").write("\n".join(lines) + "\n")
        if charge:
            open(os.path.join(d, ".CHRG"), "w").write(f"{charge}\n")
        if uhf:
            open(os.path.join(d, ".UHF"), "w").write(f"{uhf}\n")
        r = subprocess.run([GXTB_V1, "-c", "coord"], cwd=d, capture_output=True, text=True,
                           timeout=600)
        recs = read_records(os.path.join(d, "gxtbrestart"))
        import re
        m = re.search(r"nsao\s+(\d+)", r.stdout)
        nsao = int(m.group(1))
        pops = [{"sym": mm[0], "nel": float(mm[3]),
                 "s": float(mm[4]), "p": float(mm[5]), "d": float(mm[6]), "f": float(mm[7])}
                for mm in re.findall(
                    r"^##\s+\d+\s+([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                    r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
                    r.stdout, re.M)]
        _f = lambda s: [float(x) for x in  # noqa: E731
                        re.findall(r"-?\d+\.\d+(?:[eEdD][+-]?\d+)?", s)]
        eps_ev = sum((_f(l) for l in re.findall(r"^ eps  :(.*)$", r.stdout, re.M)), [])
        occ = sum((_f(l) for l in re.findall(r"^ occ\. :(.*)$", r.stdout, re.M)), [])
    state = {"nsao": nsao, "records": recs, "pops": pops,
             "eps_ev": eps_ev, "occ": occ, "raw": r.stdout}
    ntri = nsao * (nsao + 1) // 2
    for rec in recs:
        vals = np.frombuffer(rec, dtype="<f8")
        if len(vals) == ntri and "P" not in state:
            state["P"] = unpack_symmetric(vals, nsao)
        elif len(vals) == nsao * nsao and "C" not in state:
            state["C"] = vals.reshape(nsao, nsao)
        elif len(vals) == 2 * ntri and "P_a" not in state:
            # UKS: the two packed spin densities back-to-back (decoded fifty-first push;
            # the empty beta channel of a 1-electron system shows as literal zeros)
            state["P_a"] = unpack_symmetric(vals[:ntri], nsao)
            state["P_b"] = unpack_symmetric(vals[ntri:], nsao)
        elif len(vals) == 2 * nsao * nsao and "C_a" not in state:
            state["C_a"] = vals[:nsao * nsao].reshape(nsao, nsao)
            state["C_b"] = vals[nsao * nsao:].reshape(nsao, nsao)
    return state


def _gate():
    import eeqbc
    import overlap
    cases = {k: eeqbc.PROBES[k] for k in ("water", "HCl", "AcCl")}
    fails = 0
    for label, (atoms, chg, uhf) in cases.items():
        st = converged_state(atoms, charge=chg, uhf=uhf)
        zs = [eeqbc._SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
        xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * BOHR
        shells, _ = overlap.build_shells(zs, xyz, charge=chg)
        S = overlap.overlap(zs, xyz, charge=chg, shells=shells, ao_order="oracle")
        if "P" not in st or S.shape[0] != st["nsao"]:
            print(f"  {label:8s} record shapes {sorted(len(r)//8 for r in st['records'])} "
                  f"nsao {st['nsao']}: no P found  MISMATCH")
            fails += 1
            continue
        PS = st["P"] @ S
        # per-shell Mulliken populations in our shell order
        off, res = 0, {}
        for sh in shells:
            w = 2 * sh["l"] + 1
            res.setdefault((sh["at"], sh["l"]), 0.0)
            res[(sh["at"], sh["l"])] += float(np.trace(PS[off:off + w, off:off + w]))
            off += w
        worst, nel = 0.0, float(np.trace(PS))
        for i, p in enumerate(st["pops"]):
            for lname, l in (("s", 0), ("p", 1), ("d", 2), ("f", 3)):
                if (i, l) in res:
                    worst = max(worst, abs(res[(i, l)] - p[lname]))
        nel_ref = sum(p["nel"] for p in st["pops"])
        ok = worst <= 2e-3 and abs(nel - nel_ref) <= 2e-3
        fails += 0 if ok else 1
        print(f"  {label:8s} Tr(PS) {nel:9.4f} (ref {nel_ref:9.4f})  "
              f"worst shell-pop diff {worst:.2e}  {'MATCH' if ok else 'MISMATCH'}")
    if fails:
        print(f"\nRESTART GATE: {fails} MISMATCH")
        sys.exit(1)
    print("\nRESTART GATE PASSED: the parsed density + OUR overlap reproduce the oracle's "
          "Mulliken populations -- P is ours to use, and the AO conventions agree")


if __name__ == "__main__":
    _gate()
