"""eeq.py -- parse the PUBLIC EEQ(BC) charge-model file (param/gxtb/eeq).

Layout (decoded from the v1.1 asset; the selftest pins it):
  - a DATE-STAMP header line ("Sun Dec 15 12:18:58 CET 2024" -- the same stamp the oracle prints
    over its EEQ block, so the binary and this file are the same vintage; it also splits into six
    tokens and will happily masquerade as a six-float parameter row if parsed naively -- measured).
  - 103 rows of 10 floats: one row per element Z = 1..103, eight meaningful parameters + two
    zero-padded trailing columns (every element row in this release ends 0.0 0.0).
  - NO global-parameter row in this file: whatever globals EEQ(BC) needs live in the main
    parameter file (../param/gxtb/gxtb_parameters) or the code.

Which column is chi/eta/CN-coupling/radius etc. is the SI's EEQ(BC) section; the mapping is
fixed the moment our implementation reproduces the oracle's printed per-atom (CN, q_CN,
CN(basis), q) block -- oracle.py already parses that block, so the check harness exists.

    from eeq import parse
    P = parse()   # {Z: [p1..p8]} for Z = 1..103
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "param", "gxtb", "eeq")


def _numeric(line):
    try:
        return [float(x) for x in line.split()]
    except ValueError:
        return None


def parse(path=DEFAULT):
    rows = [v for l in open(path) if l.strip() and (v := _numeric(l)) is not None]
    return {z: r[:8] for z, r in enumerate([r for r in rows if len(r) == 10], start=1)}


if __name__ == "__main__":
    P = parse()
    assert len(P) == 103, len(P)
    rows = [v for l in open(DEFAULT) if l.strip() and (v := _numeric(l)) is not None]
    assert all(r[8] == 0.0 and r[9] == 0.0 for r in rows if len(r) == 10)
    assert not [r for r in rows if len(r) != 10], "only 10-wide numeric rows expected"
    print(f"elements: {len(P)} (Z 1..103); header is a date stamp, no global row in this file")
    print(f"H row:  {P[1]}")
    print(f"O row:  {P[8]}")
    print("EEQ PARSER SELFTEST PASSED")
