"""fit_basis_cn.py -- decode the q-vSZP BASIS coordination number from oracle measurements.

The basis CN (the oracle's `CN(basis)` column, feeding q_eff via Eq. 28) is NOT the EEQ-BC CN
(different values) and NOT CEH's CN (tested: CEH radii give 0.86 for water's O-H pair count where
the oracle says 0.41). Its constants live somewhere in param/gxtb/gxtb_parameters. Strategy, the
increments trick again: MEASURE the model from the oracle, FIT its constants, then FIND those
floats in the parameter file to pin the layout.

Model family (mctc/SI convention): pair count = 0.5*(1 + erf(-k*(r - rc)/rc^ne)), CN = sum.
  rc from element radii, convention SUM (rc = Ri+Rj) or MEAN (rc = (Ri+Rj)/2); ne in {1, 0.75}.

Homonuclear curves (H2, N2, O2 at several distances) give exactly-determined fits per element;
heteronuclear pairs (water, CO, HCl) then DISCRIMINATE the convention. The EEQ block prints
before the SCF, so every geometry yields data regardless of SCF behavior.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oracle  # noqa: E402


def erfinv(y, lo=-6.0, hi=6.0):
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if math.erf(mid) < y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cn_basis(atoms, charge=0, uhf=0):
    r = oracle.run(atoms, charge=charge, uhf=uhf)
    return [row["cn_basis"] for row in r["eeq"]]


def main():
    BOHR = 1.8897261254578281
    # homonuclear scans (distances in Angstrom for oracle.run)
    scans = {
        "H": [1.2, 1.4, 1.6, 2.0, 2.6],
        "N": [1.7, 2.074, 2.4, 3.0],
        "O": [1.9, 2.28, 2.7, 3.2],
    }
    uhfs = {"H": 0, "N": 0, "O": 2}
    data = {}                       # sym -> [(r_bohr, pair_count)]
    for sym, dists in scans.items():
        rows = []
        for d_bohr in dists:
            d_ang = d_bohr / BOHR
            c = cn_basis([(sym, 0.0, 0.0, 0.0), (sym, 0.0, 0.0, d_ang)], uhf=uhfs[sym])
            rows.append((d_bohr, c[0]))
            print(f"  {sym}2 @ {d_bohr:5.3f} bohr  pair count {c[0]:.4f}")
        data[sym] = rows

    # fit per convention: y = erfinv(2c-1) = -k*(r-rc)/rc^ne  ->  linear in r once rc known;
    # scan rc per element on a grid, k must come out GLOBAL (same for all elements) -- that is
    # the model's own consistency check.
    best = None
    for ne in (1.0, 0.75):
        # for each element, fit (k_elem, rc_elem) by nested 1-D search on rc
        fits = {}
        for sym, rows in data.items():
            best_e = None
            rc0 = 0.05
            while rc0 < 6.0:
                ys = [erfinv(max(-0.999999, min(0.999999, 2 * c - 1))) for _r, c in rows]
                xs = [-(r - rc0) / rc0**ne for r, _c in rows]
                sxx = sum(x * x for x in xs)
                k = sum(x * y for x, y in zip(xs, ys)) / sxx if sxx else 0.0
                res = sum((y - k * x) ** 2 for x, y in zip(xs, ys))
                if best_e is None or res < best_e[2]:
                    best_e = (rc0, k, res)
                rc0 += 0.0005
            fits[sym] = best_e
        ks = [k for _rc, k, _res in fits.values()]
        spread = max(ks) - min(ks)
        tot_res = sum(res for _rc, _k, res in fits.values())
        print(f"  ne={ne}: per-element k = " +
              ", ".join(f"{s}:{k:.4f}" for s, (_rc, k, _r) in fits.items()) +
              f"  (spread {spread:.4f}, res {tot_res:.2e})")
        if best is None or tot_res < best[2]:
            best = (ne, fits, tot_res)

    ne, fits, _ = best
    kglob = sum(k for _rc, k, _res in fits.values()) / len(fits)
    print(f"\n  BEST: ne={ne}, global k ~ {kglob:.4f}")
    for sym, (rc, k, res) in fits.items():
        print(f"    {sym}: rc(homo pair) = {rc:.4f} bohr  -> R({sym}) = {rc/2:.4f} (sum conv)"
              f" or {rc:.4f} (mean conv)")

    # heteronuclear discrimination: water O-H pair count = 0.41209 at r=1.80885
    r_oh = 1.80885
    c_target = 0.41209
    y = erfinv(2 * c_target - 1)
    for name, rc in (("SUM  R_O/2+R_H/2 -> (rcO+rcH)/2", (fits["O"][0] + fits["H"][0]) / 2),):
        pred = 0.5 * (1 + math.erf(-kglob * (r_oh - rc) / rc**ne))
        print(f"  water O-H check [{name}]: rc={rc:.4f} -> count {pred:.4f} (oracle {c_target})")

    # grep the fitted radii in the main parameter file
    import re
    text = open(os.path.join(HERE, "..", "param", "gxtb", "gxtb_parameters")).read()
    print("\n  searching gxtb_parameters for fitted radii (leading digits):")
    for sym, (rc, _k, _res) in fits.items():
        for cand, tag in ((rc / 2, "R=rc/2"), (rc, "R=rc")):
            frag = f"{cand:.3f}"[:5]
            hits = [m.start() for m in re.finditer(re.escape(frag), text)]
            print(f"    {sym} {tag} {cand:.4f}: '{frag}' -> {len(hits)} hit(s)")


if __name__ == "__main__":
    main()
