"""The U(q)-fold test on the NEUTRALS (no new oracle runs). Hypothesis: the binary's
second-order kernels run at charge-updated Hubbards U_l(q) = U_l + c*Gamma_A*q_A (chain
variants below), so the Euler-FD Gamma-part on a NEUTRAL (where the Q-weighted explicit
offsite is zero) is entirely dES2/dGamma = c*q_chain*sum dES2/dU. The decoded ES2
kernels (onsite s*harm, offsite plain KO) give dES2/dU analytically at the measured
Mulliken charges. One constant c, four measurements (OH. and H2O, sides H and O):
solve c on OH. H-side, PREDICT the other three. Variants: chain q = qat vs q_l;
Gamma bare vs k3Gs-dressed."""
import math
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402

BOHR = K.BOHR
GAM = {1: 0.8141717488, 8: 0.1454}
U = {1: [1.0062], 8: [1.4936, 0.8473]}
NSH = {1: 1, 8: 2}
SR = {1: 0.4726, 8: 0.08247 * 6 + 0.0920}
SYM = {1: "H", 8: "O"}
K3GS = -0.05708

MEAS = {("OH.", 1): -0.002664, ("OH.", 8): -0.000026,
        ("H2O", 1): -0.004420, ("H2O", 8): -0.000876}

ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
W = [[0.0, 0.0, 0.0],
     [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
     [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
CASES = [("OH.", [8, 1], [[0, 0, 0], [0, 0, 1.83]], 0, 1),
         ("H2O", [8, 1, 1], W, 0, 0)]


def shell_q(zs, xyz, charge, uhf):
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    rec = fock_recon.fock_ao(atoms, charge=charge, uhf=uhf)
    st = rec["state"]
    P = st["P"] if "P" in st else st["P_a"] + st["P_b"]
    S = rec["S"]
    net = np.diag(P @ S)
    meta = []
    for at, z in enumerate(zs):
        meta.append((at, 0))
        if NSH[z] > 1:
            meta += [(at, 1)] * 3
    q = {}
    for nn, (at, l) in zip(net, meta):
        q[(at, l)] = q.get((at, l), 0.0) + nn
    for (at, l) in list(q):
        q[(at, l)] = K.REFOCC[zs[at]].get(l, 0.0) - q[(at, l)]
    return q


def dES2_dU(q, zs, Rab, at_sel, l_sel, scope="both"):
    """analytic d(ES2)/d(U_{l_sel} of atom at_sel) with the decoded kernels.
    scope: 'onsite' | 'offsite' | 'both' -- which kernels carry U(q)."""
    z = zs[at_sel]
    Us = U[z][l_sel]
    d = 0.0
    # onsite: ES2_on = (1/2) s_z sum_{l,l'} q_l q_l' harm(U_l, U_l')
    for la in range(NSH[z]) if scope in ("onsite", "both") else []:
        for lb in range(NSH[z]):
            qa, qb = q[(at_sel, la)], q[(at_sel, lb)]
            Ua, Ub = U[z][la], U[z][lb]
            dh = 0.0
            if la == l_sel:
                dh += 2 * Ub * Ub / (Ua + Ub) ** 2
            if lb == l_sel:
                dh += 2 * Ua * Ua / (Ua + Ub) ** 2
            if la == l_sel and lb == l_sel:
                dh = 1.0        # harm(U,U) = U
            d += 0.5 * SR[z] * qa * qb * dh
    # offsite: ES2_off = sum_{B != A} sum_{l'} q_lA q_l'B gammaKO; dgamma/dU =
    # gamma^2/(2U^2)
    for (bt, lb), qb in (q.items() if scope in ("offsite", "both") else []):
        if bt == at_sel:
            continue
        R = Rab[at_sel][bt]
        g = 1.0 / (R + 0.5 * (1.0 / Us + 1.0 / U[zs[bt]][lb]))
        d += q[(at_sel, l_sel)] * qb * g * g / (2 * Us * Us)
    return d


print("U(q)-fold sensitivities S_X = sum_l q_chain * dES2/dU_l  (per variant):")
rows = {}
rows_scoped = {}
for name, zs, xyz, chg, uhf in CASES:
    q = shell_q(zs, xyz, chg, uhf)
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    nat = len(zs)
    Rab = [[float(np.linalg.norm(np.array(xyz[a]) - np.array(xyz[b])))
            for b in range(nat)] for a in range(nat)]
    for z_sel in sorted(set(zs)):
        s_at = 0.0      # chain q = qat
        s_l = 0.0       # chain q = q_l
        s_sc = {"onsite": [0.0, 0.0], "offsite": [0.0, 0.0]}
        for at in range(nat):
            if zs[at] != z_sel:
                continue
            for l in range(NSH[z_sel]):
                dd = dES2_dU(q, zs, Rab, at, l)
                s_at += qat[at] * dd
                s_l += q[(at, l)] * dd
                for sc in ("onsite", "offsite"):
                    ds = dES2_dU(q, zs, Rab, at, l, sc)
                    s_sc[sc][0] += qat[at] * ds
                    s_sc[sc][1] += q[(at, l)] * ds
        rows[(name, z_sel)] = (s_at, s_l)
        rows_scoped[(name, z_sel)] = s_sc
        print(f"  {name} {SYM[z_sel]}: S(qat-chain) {s_at:+.6f}  "
              f"S(ql-chain) {s_l:+.6f}   measured extra {MEAS[(name, z_sel)]:+.6f}")

print("\nsolve c on OH. H-side, PREDICT the rest:")
for chain, idx in (("qat", 0), ("ql", 1)):
    for gtag, gfac in (("bare-G", 1.0), ("k3Gs-dressed", K3GS)):
        s0 = rows[("OH.", 1)][idx]
        c = MEAS[("OH.", 1)] / (GAM[1] * gfac * s0)
        print(f"  chain={chain:3s} {gtag:12s}: c = {c:+.5f}")
        ok = True
        for key in (("OH.", 8), ("H2O", 1), ("H2O", 8)):
            s = rows[key][idx]
            pred = c * GAM[key[1]] * gfac * s
            meas = MEAS[key]
            hit = abs(pred - meas) <= max(0.001, 0.25 * abs(meas))
            ok = ok and hit
            print(f"      {key[0]} {SYM[key[1]]}: pred {pred:+.6f}  meas {meas:+.6f}  "
                  f"{'HIT' if hit else 'miss'}")
        print(f"      => {'FOLD CONFIRMED (this variant)' if ok else 'variant fails'}")

print()
print("SCOPED folds -- solve c on OH. H-side, predict the rest; plus FIXED c=1:")
for scope in ("onsite", "offsite"):
    for chain, idx in (("qat", 0), ("ql", 1)):
        s0 = rows_scoped[("OH.", 1)][scope][idx]
        if abs(s0) < 1e-12:
            continue
        c = MEAS[("OH.", 1)] / (GAM[1] * s0)
        line = [f"scope={scope:7s} chain={chain:3s} bare-G: c = {c:+.5f}"]
        ok = True
        for key in (("OH.", 8), ("H2O", 1), ("H2O", 8)):
            sv = rows_scoped[key][scope][idx]
            pred = c * GAM[key[1]] * sv
            meas = MEAS[key]
            hit = abs(pred - meas) <= max(0.001, 0.25 * abs(meas))
            ok = ok and hit
            line.append(f"    {key[0]} {SYM[key[1]]}: pred {pred:+.6f} meas "
                        f"{meas:+.6f} {'HIT' if hit else 'miss'}")
        line.append(f"    => {'FOLD CONFIRMED' if ok else 'variant fails'}")
        print(chr(10).join(line))
print()
print("FIXED c = +1 (Gamma IS dU/dq), offsite-only, bare-G -- pure prediction:")
for key in (("OH.", 1), ("OH.", 8), ("H2O", 1), ("H2O", 8)):
    for chain, idx in (("qat", 0), ("ql", 1)):
        sv = rows_scoped[key][scope := "offsite"][idx]
        pred = 1.0 * GAM[key[1]] * sv
        print(f"  {key[0]} {SYM[key[1]]} ({chain}): pred {pred:+.6f}  "
              f"meas {MEAS[key]:+.6f}")
