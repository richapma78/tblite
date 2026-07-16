"""oracle.py -- drive the v1.1 g-xTB binary and parse EVERY intermediate it prints.

The v1 oracle (public binary + public parameter files, see ../GXTB-PLAN.md) is this project's
ground truth. It is unusually talkative, and that is the whole G1 strategy: it prints the EEQ(BC)
charges and coordination numbers, the q-vSZP per-atom adaptation (with its functional form in the
header), the basis dimensions, the SCF trace, eigenvalues, shell populations, Wiberg bond orders,
Janak IP/EA, and the energy DECOMPOSED BY TERM (ES1 / ES2+3 / multipole / Mulliken exchange /
spin / electronic / atomic increments / dispersion / repulsion). Every one of those is a
checkpoint our implementation must hit BEFORE the total is allowed to match.

    from oracle import run
    r = run([("O", 0.0, 0.0, 0.117), ("H", 0.0, 0.757, -0.471), ("H", 0.0, -0.757, -0.471)])
    r["terms"]["electronic"], r["eeq"][0]["q"], r["gap_ev"], r["total"]

Numerical gradients: run(..., grad=True) adds r["gradient"] (Eh/Bohr, per atom) parsed from the
Turbomole-style gradient file the -grad flag writes. Slow (finite differences) but they are the
reference our future ANALYTIC gradients (SI eqs.; gate G4) get checked against.
"""
import os
import re
import subprocess
import tempfile

GXTB_V1 = os.environ.get("GXTB_V1_EXE", "/opt/gxtb-v1/gxtb")
BOHR = 1.8897261254578281

_TERMS = ("ES1 (charge SIE)", "ES2+3", "ES multipole", "ES total", "Ex (Mulliken)", "Espinpol",
          "electronic", "atomic core increments", "dispersion", "nuclear repulsion")


def _floats(s):
    return [float(x) for x in re.findall(r"-?\d+\.\d+(?:[eEdD][+-]?\d+)?", s)]


def run(atoms, charge=0, uhf=0, grad=False, timeout=600):
    """atoms: [(symbol, x, y, z), ...] in ANGSTROM. Returns the parsed intermediate stack."""
    lines = ["$coord"]
    for sym, x, y, z in atoms:
        lines.append(f" {x * BOHR:.10f} {y * BOHR:.10f} {z * BOHR:.10f} {sym.lower()}")
    lines.append("$end")
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "coord"), "w").write("\n".join(lines) + "\n")
        if charge:
            open(os.path.join(d, ".CHRG"), "w").write(f"{charge}\n")
        if uhf:
            open(os.path.join(d, ".UHF"), "w").write(f"{uhf}\n")
        cmd = [GXTB_V1, "-c", "coord"] + (["-grad"] if grad else [])
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=timeout)
        out = r.stdout
        res = {"raw": out, "atoms": atoms, "charge": charge, "uhf": uhf}

        # EEQ(BC) block: "    1 O    1.1524  -0.3590   0.8242  -0.2803"
        res["eeq"] = [{"sym": m[1], "cn": float(m[2]), "q_cn": float(m[3]),
                       "cn_basis": float(m[4]), "q": float(m[5])}
                      for m in re.findall(
                          r"^\s+(\d+)\s+([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                          r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", out, re.M)]
        # q-vSZP AO adaptation: "   8    -0.28031     0.82418    -0.28124  0.28480 -0.04148 -0.03793"
        res["ao_adapt"] = [{"z": int(m[0]), "q": float(m[1]), "cn": float(m[2]),
                            "dq": float(m[3]), "dcn": float(m[4]), "dqcn": float(m[5]),
                            "total": float(m[6])}
                           for m in re.findall(
                               r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                               r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", out, re.M)]
        for key, pat in (("ncao", r"ncao\s+(\d+)"), ("nsao", r"nsao\s+(\d+)"),
                         ("npr", r"npr\s+(\d+)"), ("acpsao", r"acpsao\s+(\d+)"),
                         ("nel", r"nel\s+(\d+)"), ("nalpha", r"nalpha\s+(\d+)"),
                         ("nbeta", r"nbeta\s+(\d+)")):
            m = re.search(pat, out)
            res[key] = int(m.group(1)) if m else None
        res["uks"] = bool(re.search(r"UKS \?\s+T", out))
        m = re.search(r"(\d+)\s+scf iterations", out)
        res["scf_iterations"] = int(m.group(1)) if m else None
        m = re.search(r"gap \(eV\)\s*:\s*(-?\d+\.\d+)", out)
        res["gap_ev"] = float(m.group(1)) if m else None
        m = re.search(r"Janaks theorem for IP\s*:\s*(-?\d+\.\d+)", out)
        res["janak_ip_ev"] = float(m.group(1)) if m else None
        m = re.search(r"Janaks theorem for EA\s*:\s*(-?\d+\.\d+)", out)
        res["janak_ea_ev"] = float(m.group(1)) if m else None
        # eigenvalues: "eps  :   -34.8179  -13.8363 ..." lines follow "occ. :" lines
        res["occ"] = sum((_floats(l) for l in re.findall(r"^ occ\. :(.*)$", out, re.M)), [])
        res["eps_ev"] = sum((_floats(l) for l in re.findall(r"^ eps  :(.*)$", out, re.M)), [])
        # shell populations: "##     1 O   1.99  6.0  6.6558   1.784  4.871  0.000  0.000"
        res["pops"] = [{"sym": m[0], "cn": float(m[1]), "z": float(m[2]), "nel": float(m[3]),
                        "s": float(m[4]), "p": float(m[5]), "d": float(m[6]), "f": float(m[7])}
                       for m in re.findall(
                           r"^##\s+\d+\s+([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                           r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                           r"\s+(-?\d+\.\d+)", out, re.M)]
        res["terms"] = {}
        for t in _TERMS:
            m = re.search(rf"^\s*{re.escape(t)}\s*:\s*(-?\d+\.\d+)", out, re.M)
            if m:
                res["terms"][t] = float(m.group(1))
        m = re.search(r"^\s*total\s+(-?\d+\.\d+)", out, re.M)
        res["total"] = float(m.group(1)) if m else None
        if grad:
            gpath = os.path.join(d, "gradient")
            res["gradient"] = None
            if os.path.exists(gpath):
                gl = open(gpath).read()
                rows = re.findall(r"^\s*(-?\d+\.\d+[dDeE]?[+-]?\d*)\s+(-?\d+\.\d+[dDeE]?[+-]?\d*)"
                                  r"\s+(-?\d+\.\d+[dDeE]?[+-]?\d*)\s*$", gl, re.M)
                if rows:
                    res["gradient"] = [[float(x.replace("D", "E").replace("d", "e"))
                                        for x in row] for row in rows[-len(atoms):]]
        return res


if __name__ == "__main__":
    r = run([("O", 0.0, 0.0, 0.117), ("H", 0.0, 0.757, -0.471), ("H", 0.0, -0.757, -0.471)])
    assert r["total"] is not None and abs(r["total"] - (-76.43743981)) < 5e-6, r["total"]
    assert len(r["eeq"]) == 3 and abs(r["eeq"][0]["q"] - (-0.2803)) < 2e-3, r["eeq"]
    assert r["npr"] == 40 and r["nel"] == 8 and len(r["terms"]) >= 9, (r["npr"], r["terms"])
    print(f"water: total {r['total']}, {len(r['eeq'])} EEQ rows, npr {r['npr']}, "
          f"{len(r['terms'])} terms, gap {r['gap_ev']:.2f} eV")
    print("ORACLE DRIVER SELFTEST PASSED")
