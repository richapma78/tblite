"""extract_tables.py -- pull the data tables EEQ(BC) needs out of the Fortran sources, once.

Sources (sibling clones, exact revisions our reproduction targets):
  - multicharge/src/multicharge/param/eeqbc2025.f90 : the PUBLISHED eeqbc2025 parameter arrays
    (chi, eta, rad, kcnchi, kqchi, kqeta, cap, cov_radii, avg_cn). Our fork's param/gxtb/eeq file
    carries 8 of these -- this script CROSS-CHECKS all 103 elements of every shared column, so a
    single transcription slip anywhere would fail loudly.
  - mctc-lib/src/mctc/data/vdwrad.f90 : pairwise van-der-Waals radii, packed lower triangle,
    literals in Angstrom (the Fortran declares aatoau*[...]).
  - mctc-lib/src/mctc/data/paulingen.f90 : Pauling electronegativities (Z<=118).

Writes prototype/data/mctc-tables.json. Run from anywhere; paths derive from this file.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MC = os.environ.get("MULTICHARGE_SRC", "C:/Projects/multicharge/src/multicharge")
MCTC = os.environ.get("MCTC_SRC", "C:/Projects/mctc-lib/src/mctc")
if not os.path.isdir(MC):                       # WSL path fallback
    MC, MCTC = "/mnt/c/Projects/multicharge/src/multicharge", "/mnt/c/Projects/mctc-lib/src/mctc"

FLOAT = re.compile(r"(-?\d+\.\d+(?:[eE][+-]?\d+)?)_wp")


def fortran_array(path, name, count):
    """All _wp literals between `:: name(...) = [` and the closing `]`."""
    text = open(path).read()
    m = re.search(rf"::\s*{name}\(.*?\)\s*=\s*(?:[\w.]+\s*\*\s*)?\[(.*?)\]", text, re.S)
    if not m:
        raise KeyError(f"{name} not found in {path}")
    vals = [float(x) for x in FLOAT.findall(m.group(1))]
    if len(vals) != count:
        raise ValueError(f"{name}: expected {count} values, got {len(vals)}")
    return vals


def main():
    n = 103
    ee = os.path.join(MC, "param", "eeqbc2025.f90")
    pub = {a: fortran_array(ee, f"eeqbc_{a}", n)
           for a in ("chi", "eta", "rad", "kcnchi", "kqchi", "kqeta", "cap",
                     "cov_radii", "avg_cn")}
    vdw = fortran_array(os.path.join(MCTC, "data", "vdwrad.f90"), "vdwrad", n * (n + 1) // 2)
    en = fortran_array(os.path.join(MCTC, "data", "paulingen.f90"), "pauling_en", 118)

    # cross-check the fork's eeq file against the published arrays, ALL elements, ALL 8 columns
    sys.path.insert(0, HERE)
    import eeq as eeqmod
    file_rows = eeqmod.parse()
    # NOTE: the extractor strips the Fortran declaration's 0.5_wp* prefix, so pub["cov_radii"]
    # holds the RAW literals -- exactly what the eeq file stores. The WORKING covalent radius
    # (used in the CN counting function) is 0.5 * raw, applied in eeqbc.py, in Bohr.
    colmap = ["chi", "eta", "rad", "kcnchi", "cov_radii", "kqeta", "kqchi", "cap"]
    worst = 0.0
    for z in range(1, n + 1):
        for c, name in enumerate(colmap):
            worst = max(worst, abs(file_rows[z][c] - pub[name][z - 1]))
    print(f"eeq file vs published eeqbc2025, all 103 elements x 8 columns: "
          f"max |diff| = {worst:.2e}  {'IDENTICAL' if worst < 1e-9 else 'MISMATCH'}")
    assert worst < 1e-9, "the v1 eeq file is NOT the published eeqbc2025 set -- stop and look"

    out = {"_what": "data tables for the EEQ(BC) reproduction (see extract_tables.py header)",
           "_units": {"vdwrad_pair": "Angstrom literals (Fortran multiplies aatoau; "
                                     "multicharge then multiplies autoaa back)",
                      "cov_radii": "published array is 0.5*raw, in Bohr",
                      "pauling_en": "raw Pauling; eeqbc2025 normalizes /3.98 + actinide patches"},
           "eeqbc2025": pub, "vdwrad_packed_lower_angstrom": vdw, "pauling_en": en}
    path = os.path.join(HERE, "data", "mctc-tables.json")
    json.dump(out, open(path, "w"))
    print(f"wrote {path} ({os.path.getsize(path)//1024} KB)")


if __name__ == "__main__":
    main()
