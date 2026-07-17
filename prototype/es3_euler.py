"""ES3 offsite via the Euler-FD instrument. E3 is linear in each Gamma_X, so
Gamma_X * dES23/dGamma_X = the X-part of E3 EXACTLY (central FD; no ES2 model involved).
Subtracting the decoded onsite (tau=-1/2U^2, k3Gs=-0.05708, k3Gp=-0.0668) leaves CLEAN
per-element offsite targets. Solve (k3, k3x) for bare-Gamma variants 1/2 on the H3+ and
HeH+ ladders; OH- is held out as the out-of-sample check.
PRE-DECLARED: the offsite is CLOSED if solve rms <= 0.003 with no point over 0.006 AND
the OH- out-of-sample parts land within 0.008; else it stays open and says so."""
import math
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import oracle  # noqa: E402
import params  # noqa: E402

BOHR = K.BOHR
K3GS, K3GP = -0.05708, -0.06679
GAM = {1: 0.8141717488, 2: 5.0001493881, 8: 0.1454, 9: 0.567856801}
U = {1: [1.0062], 2: [0.6537], 8: [1.4936, 0.8473], 9: [1.5589, 0.7561]}
NSH = {1: 1, 2: 1, 8: 2, 9: 2}
SYM = {1: "H", 2: "He", 8: "O", 9: "F"}


def tri(R):
    return [[0, 0, 0], [R, 0, 0], [R / 2, R * math.sqrt(3) / 2, 0]]


def es23(atoms, charge):
    r = oracle.run(atoms, charge=charge)
    return float(re.search(r"^\s*ES2\+3\s*:\s*(-?\d+\.\d+)", r["raw"], re.M).group(1))


def shell_q(zs, xyz, charge):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge)
    P, S = rec["state"]["P"], rec["S"]
    m = np.diag((P / 2.0) @ S)
    meta = []
    for at, z in enumerate(zs):
        meta.append((at, 0))
        if NSH[z] > 1:
            meta += [(at, 1)] * 3
    q = {}
    for mm, (at, l) in zip(m, meta):
        q[(at, l)] = q.get((at, l), 0.0) + 2 * mm
    for (at, l) in list(q):
        q[(at, l)] = K.REFOCC[zs[at]].get(l, 0.0) - q[(at, l)]
    return q


def e3_onsite_part(q, zs, at_sel):
    """onsite E3 terms living on atom at_sel (all carry Gamma of that element)."""
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    e3 = 0.0
    for (at, la), qa in q.items():
        if at != at_sel:
            continue
        z = zs[at]
        for lb in range(NSH[z]):
            qb = q[(at, lb)]
            ta = -1.0 / (2 * U[z][la] ** 2)
            tb = -1.0 / (2 * U[z][lb] ** 2)
            ga = (K3GS if la == 0 else K3GP) * GAM[z]
            gb = (K3GS if lb == 0 else K3GP) * GAM[z]
            e3 += (1 / 6) * qa * qb * qat[at] * (ta * ga + tb * gb)
    return e3


def e3_off_part(q, zs, Rab, k3, k3x, var, z_sel):
    """offsite E3 terms carrying Gamma of element z_sel (the qat*tau*G(at) halves
    where zs[at] == z_sel)."""
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    e3 = 0.0
    for (at, la), qa in q.items():
        if zs[at] != z_sel:
            continue
        for (bt, lb), qb in q.items():
            if at == bt:
                continue
            R = Rab[at][bt]
            ub = 0.5 * (U[zs[at]][la] + U[zs[bt]][lb])
            w = ub * ub if var == 1 else ub
            ex = math.exp(-k3x * w * R)
            tau = k3 * ub * R * ex * (1 - k3x * w * R)
            e3 += (1 / 6) * qa * qb * qat[at] * tau * GAM[zs[at]]
    return e3


CASES = [("H3+@1.5", [1, 1, 1], tri(1.5), 1),
         ("H3+@1.65", [1, 1, 1], tri(1.65), 1),
         ("H3+@1.9", [1, 1, 1], tri(1.9), 1),
         ("HeH+@1.46", [2, 1], [[0, 0, 0], [0, 0, 1.46]], 1),
         ("HeH+@2.0", [2, 1], [[0, 0, 0], [0, 0, 2.0]], 1),
         ("HeH+@2.5", [2, 1], [[0, 0, 0], [0, 0, 2.5]], 1),
         ("OH-@1.83", [8, 1], [[0, 0, 0], [0, 0, 1.83]], -1)]

DELTA = 0.05
obs = []          # (name, z_sel, q, zs, Rab, clean_offsite_target, holdout)
P0 = params.parse()
try:
    for name, zs, xyz, chg in CASES:
        atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
                 for z, (x, y_, zc) in zip(zs, xyz)]
        q = shell_q(zs, xyz, chg)
        nat = len(zs)
        Rab = [[float(np.linalg.norm(np.array(xyz[a]) - np.array(xyz[b])))
                for b in range(nat)] for a in range(nat)]
        for z_sel in sorted(set(zs)):
            g0 = GAM[z_sel]
            params.write_perturbed({(z_sel, 1, 6): g0 + DELTA})
            ep = es23(atoms, chg)
            params.write_perturbed({(z_sel, 1, 6): g0 - DELTA})
            em = es23(atoms, chg)
            params.write_perturbed({})
            part = g0 * (ep - em) / (2 * DELTA)
            onsite = sum(e3_onsite_part(q, zs, at)
                         for at in range(nat) if zs[at] == z_sel)
            off = part - onsite
            holdout = name.startswith("OH-")
            obs.append((name, z_sel, q, zs, Rab, off, holdout))
            print(f"  {name} G_{SYM[z_sel]}: E3-part {part:+.5f}  onsite {onsite:+.5f}"
                  f"  CLEAN offsite {off:+.5f}{'   [holdout]' if holdout else ''}",
                  flush=True)
finally:
    params.write_perturbed({})

fit = [o for o in obs if not o[6]]
for var in (1, 2):
    best = None
    for k3 in np.arange(-8.0, 8.01, 0.2):
        if abs(k3) < 1e-9:
            continue
        for k3x in np.arange(0.02, 3.01, 0.04):
            r2 = sum((e3_off_part(q, zs, Rab, k3, k3x, var, z) - off) ** 2
                     for _, z, q, zs, Rab, off, _h in fit)
            if best is None or r2 < best[0]:
                best = (r2, k3, k3x)
    r2, k3, k3x = best
    for _ in range(4):
        for k3t in np.arange(k3 - 0.2, k3 + 0.2, 0.02):
            for k3xt in np.arange(max(0.01, k3x - 0.06), k3x + 0.06, 0.005):
                r2t = sum((e3_off_part(q, zs, Rab, k3t, k3xt, var, z) - off) ** 2
                          for _, z, q, zs, Rab, off, _h in fit)
                if r2t < r2:
                    r2, k3, k3x = r2t, k3t, k3xt
    rms = math.sqrt(r2 / len(fit))
    worst = max(abs(e3_off_part(q, zs, Rab, k3, k3x, var, z) - off)
                for _, z, q, zs, Rab, off, _h in fit)
    print(f"\nvariant {var}: k3 = {k3:+.4f}  k3x = {k3x:.4f}  rms {rms:.5f}  "
          f"worst {worst:.5f}")
    for name, z, q, zs, Rab, off, h in obs:
        pred = e3_off_part(q, zs, Rab, k3, k3x, var, z)
        tag = "  [OUT-OF-SAMPLE]" if h else ""
        print(f"    {name} G_{SYM[z]}: clean {off:+.5f}  pred {pred:+.5f}  "
              f"d {pred - off:+.5f}{tag}")
    oos_ok = all(abs(e3_off_part(q, zs, Rab, k3, k3x, var, z) - off) <= 0.008
                 for _, z, q, zs, Rab, off, h in obs if h)
    closed = rms <= 0.003 and worst <= 0.006 and oos_ok
    print(f"  PRE-DECLARED VERDICT: {'CLOSED' if closed else 'OPEN'}")
