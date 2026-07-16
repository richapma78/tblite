"""repulsion.py -- Term 5: the semi-classical charge-dependent repulsion (SI Sec. 1.6), gated.

  E_rep = 1/2 sum_{A!=B} Zeff_A * Zeff_B * f(R_AB)
  Zeff_A = zeff0_A * (1 - kq_A*q_A - kq2_A*q_A^2)   [BOTH signs measured by FD probing: ratio -1.000 exposed the + reading as wrong]          q = the EEQ(BC) charge (Term 2)
  f(R)   = (1 + p1_AB/(2R) + kpen2/R^2 + kpen3/R^3 + kpen4/R^4) * exp(-a_AB*(R + s*R0_AB)^1.5)
  p1_AB  = kpen1_A + kpen1_B  (class value: H/He get kpen1_hhe, others kpen1)
  a_A    = alpha0_A * (1 + kcn_rep_A * sqrt(CN_A))          CN = Eq. 47 (internal CN)
  a_AB   = a_A*a_B/(a_A+a_B) ; R0_AB = sqrt(r0_A*r0_B)      (SI combination rules)
  CN_A   = sum_B 0.5*(1 + erf(kcn_glob*(R-rc)/rc)), rc = rcov_cn_A + rcov_cn_B  [sign convention
           under test: kcn_glob is stored NEGATIVE; rc as sum vs arithmetic mean is Eq.47's
           "arithmetical average" -- BOTH variants are gated below and the winner recorded]

The exponent OFFSET SIGN (R + R0 vs R - R0) and the rc convention are the two readings the
garbled SI text leaves open; the gate decides them empirically. Everything else is pinned by
the perturbation map (params.py docstring).

GATE (declared): reproduce the oracle's printed "nuclear repulsion" to <= 5e-7 Eh on the H2
curve (10 distances), HCl, NH4+ (charged: kq/kq2 live), water, AcCl, PdCl2 -- one variant
combination must pass ALL; if none does, record the failure, do not tune constants.
"""
import itertools
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import params  # noqa: E402

BOHR = 1.8897261254578281


def cn_eq47(zs, xyz, P, mean=False):
    k = P["globals"]["kcn_glob"]
    n = len(zs)
    cn = np.zeros(n)
    for i in range(n):
        for j in range(i):
            ri = P["element"][zs[i]]["rcov_cn"]
            rj = P["element"][zs[j]]["rcov_cn"]
            rc = 0.5 * (ri + rj) if mean else (ri + rj)
            r = float(np.linalg.norm(xyz[i] - xyz[j]))
            c = 0.5 * (1.0 + math.erf(k * (r - rc) / rc))
            cn[i] += c
            cn[j] += c
    return cn


def energy(zs, xyz, q, P=None, sign=+1, mean_rc=False, comb="harmonic"):
    """E_rep in Eh. zs: atomic numbers; xyz: (n,3) Bohr; q: EEQ(BC) charges.
    comb: alpha_AB combination -- 'gauss' = aA*aB/(aA+aB) (SI Eq. 55 as extracted),
    'harmonic' = 2*aA*aB/(aA+aB), 'geometric' = sqrt(aA*aB). The H2 tail refutes 'gauss'
    empirically (decays 150x too slowly); the extraction likely lost a factor of 2."""
    if P is None:
        P = params.parse()
    G = P["globals"]
    cn = cn_eq47(zs, xyz, P, mean=mean_rc)
    n = len(zs)
    zeff, aa, r0, p1 = np.empty(n), np.empty(n), np.empty(n), np.empty(n)
    for i, z in enumerate(zs):
        e = P["element"][z]
        zeff[i] = e["zeff0"] * (1.0 - e["kq_rep"] * q[i] - e["kq2_rep"] * q[i] ** 2)
        aa[i] = e["alpha0"] * (1.0 + e["kcn_rep"] * math.sqrt(max(cn[i], 0.0)))
        r0[i] = e["r0"]
        p1[i] = G["kpen1_hhe"] if z in (1, 2) else G["kpen1"]
    E = 0.0
    for i in range(n):
        for j in range(i):
            r = float(np.linalg.norm(xyz[i] - xyz[j]))
            if comb == "gauss":
                a_ab = aa[i] * aa[j] / (aa[i] + aa[j])
            elif comb == "harmonic":
                a_ab = 2.0 * aa[i] * aa[j] / (aa[i] + aa[j])
            else:
                a_ab = math.sqrt(aa[i] * aa[j])
            r0_ab = math.sqrt(r0[i] * r0[j])
            pen = (1.0 + (p1[i] + p1[j]) / (2.0 * r) + G["kpen2"] / r ** 2
                   + G["kpen3"] / r ** 3 + G["kpen4"] / r ** 4)
            E += zeff[i] * zeff[j] * pen * math.exp(-a_ab * (r + sign * r0_ab) ** params.KEXP_REP)
    return E


def _gate():
    import eeqbc
    import oracle
    P = params.parse()
    cases = []
    for rb in (1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0):
        cases.append((f"H2@{rb}", [("H", 0, 0, 0), ("H", 0, 0, rb / BOHR)], 0, 0))
    for name in ("HCl", "NH4+", "water", "AcCl", "PdCl2"):
        atoms, chg, uhf = eeqbc.PROBES[name]
        cases.append((name, atoms, chg, uhf))

    refs, geoms = [], []
    for name, atoms, chg, uhf in cases:
        r = oracle.run(atoms, charge=chg, uhf=uhf)
        zs = [eeqbc._SYM2Z[s.lower()] for s, _x, _y, _z in atoms]
        xyz = np.array([[x, y, z] for _s, x, y, z in atoms]) * BOHR
        q = eeqbc.charges(zs, xyz, charge=chg)["q"]
        refs.append(r["terms"]["nuclear repulsion"])
        geoms.append((name, zs, xyz, q))

    best = None
    for sign, mean_rc, comb in itertools.product((+1, -1), (False, True),
                                                 ("gauss", "harmonic", "geometric")):
        worst, wname = 0.0, ""
        for (name, zs, xyz, q), ref in zip(geoms, refs):
            e = energy(zs, xyz, q, P, sign=sign, mean_rc=mean_rc, comb=comb)
            if abs(e - ref) > worst:
                worst, wname = abs(e - ref), name
        tag = (f"exp(-a(R{'+' if sign > 0 else '-'}R0)^1.5), "
               f"rc={'mean' if mean_rc else 'sum'}, a={comb}")
        print(f"  {tag:52s} worst |dE| = {worst:.2e} ({wname})")
        if best is None or worst < best[0]:
            best = (worst, sign, mean_rc, comb)
    worst, sign, mean_rc, comb = best
    if worst > 5e-7:
        print(f"\nREPULSION GATE: FAILED (best variant worst |dE| {worst:.2e}) -- "
              "a structural reading is still wrong; do not tune.")
        sys.exit(1)
    print(f"\nREPULSION GATE PASSED: sign={'+' if sign > 0 else '-'}R0, "
          f"rc={'mean' if mean_rc else 'sum'}, a_AB={comb}; worst |dE| = {worst:.2e} over "
          f"{len(cases)} cases incl. charged NH4+ and PdCl2")


if __name__ == "__main__":
    _gate()
