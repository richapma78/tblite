"""energy_parity.py -- THE GAP LEDGER. For each gate neutral: every printed term of
the binary's energy decomposition vs ours-where-implemented, in mEh. Terms we cannot
compute yet print as MISSING (that is the point: the ranked defect list for the
substitution-grade campaign). Also measures the per-element atomic core increments
(5 atom runs) since those are constants the port must carry."""
import math
import re
import sys

import numpy as np

sys.path.insert(0, "/mnt/c/Projects/tblite-gxtb/prototype")
import constants as K  # noqa: E402
import fock_recon  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import oracle  # noqa: E402
import repulsion  # noqa: E402
import scf_h2  # noqa: E402
import params  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
P_ = params.parse()

# ---- per-element atomic core increments (ground-state atom runs)
UHF_ATOM = {1: 1, 6: 2, 7: 3, 8: 2, 9: 1}
INC = {}
print("atomic core increments (constants for the port):")
for z, u in UHF_ATOM.items():
    r = oracle.run([(SYM[z], 0, 0, 0)], uhf=u)
    m = re.search(r"^\s*atomic core increments\s*:\s*(-?\d+\.\d+)", r["raw"], re.M)
    INC[z] = float(m.group(1))
    print(f"  {SYM[z]}: {INC[z]:+.8f}")

ang = 104.5 * math.pi / 180
r_oh = 0.9572 * BOHR
W = [[0.0, 0.0, 0.0],
     [r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)],
     [-r_oh * math.sin(ang / 2), 0.0, r_oh * math.cos(ang / 2)]]
a4 = 1.087 * BOHR / math.sqrt(3)
r_nh = 1.012 * BOHR
st, ct = 0.9262, -0.3770
SYSTEMS = [
    ("H2", [1, 1], [[0, 0, 0], [0, 0, 1.4]]),
    ("F2", [9, 9], [[0, 0, 0], [0, 0, 2.668]]),
    ("HF", [1, 9], [[0, 0, 0], [0, 0, 1.733]]),
    ("H2O", [8, 1, 1], W),
    ("CH4", [6, 1, 1, 1, 1],
     [[0, 0, 0], [a4, a4, a4], [a4, -a4, -a4], [-a4, a4, -a4], [-a4, -a4, a4]]),
    ("NH3", [7, 1, 1, 1],
     [[0, 0, 0]] + [[r_nh * st * math.cos(2 * math.pi * k / 3),
                     r_nh * st * math.sin(2 * math.pi * k / 3),
                     r_nh * ct] for k in range(3)]),
]
AUFBAU = {1: {0: 1.0}, 6: {0: 2.0, 1: 2.0}, 7: {0: 2.0, 1: 3.0},
          8: {0: 2.0, 1: 4.0}, 9: {0: 2.0, 1: 5.0}}

defects = {}
for name, zs, xyz in SYSTEMS:
    atoms = [(SYM[z], x / BOHR, y_ / BOHR, zc / BOHR)
             for z, (x, y_, zc) in zip(zs, xyz)]
    try:
        rec = fock_recon.fock_ao(atoms)
        raw = rec["state"]["raw"]
    except AssertionError:
        # not fully invertible (window < nsao): printed terms only, plus the
        # density-free terms at our own SCF charges (labeled)
        rr = oracle.run(atoms)
        T = {k.strip(): float(v) for k, v in
             re.findall(r"^([A-Za-z0-9+() .]+?)\s*:\s*(-?\d+\.\d+)\s*$",
                        rr["raw"], re.M)}
        ours = {"atomic core increments": sum(INC[z] for z in zs)}
        print(f"{chr(10)}{name}: printed-vs-ours (NOT invertible; printed + "
              f"increments only)")
        for t in ["electronic", "Ex (Mulliken)", "ES1 (charge SIE)", "ES2+3",
                  "ES multipole", "ES total", "atomic core increments",
                  "dispersion", "nuclear repulsion"]:
            if t not in T:
                continue
            if t in ours:
                d = (ours[t] - T[t]) * 1000
                print(f"  {t:24s} printed {T[t]:+12.6f}  ours {ours[t]:+12.6f}  "
                      f"d {d:+9.3f} mEh")
                defects.setdefault(t, []).append(abs(d))
            else:
                print(f"  {t:24s} printed {T[t]:+12.6f}  ours     MISSING")
                defects.setdefault(t + " [MISSING]", []).append(abs(T[t]) * 1000)
        continue
    T = {k.strip(): float(v) for k, v in
         re.findall(r"^([A-Za-z0-9+() .]+?)\s*:\s*(-?\d+\.\d+)\s*$", raw, re.M)}
    P = rec["state"]["P"]
    B = GE.build(zs, np.array(xyz, float))
    S, meta, E, n = B["S"], B["meta"], B["E"], B["n"]
    mv = np.diag((P / 2.0) @ S)
    q = {}
    for i in range(n):
        at, z, l, _ = meta[i]
        q[(at, l)] = q.get((at, l), 0.0) + 2 * mv[i]
    for (at, l) in list(q):
        q[(at, l)] = E[zs[at]]["ref"][l] - q[(at, l)]
    qat = {}
    for (at, l), v in q.items():
        qat[at] = qat.get(at, 0.0) + v
    ours = {}
    # nuclear repulsion
    ours["nuclear repulsion"] = repulsion.energy(
        list(zs), np.array(xyz, float), [qat[a] for a in range(len(zs))], P_,
        sign=+1, mean_rc=True, comb="harmonic")
    # ES2+3: analytic onsite ES2 + the charge-driven block (folded KO + ES3 + Q-term)
    es2on = 0.0
    for (at, l), qa in q.items():
        z = zs[at]
        for l2 in range(E[z]["nsh"]):
            g2 = GE.SRULE[z] * 2 * E[z]["U"][l] * E[z]["U"][l2] / \
                (E[z]["U"][l] + E[z]["U"][l2])
            es2on += 0.5 * qa * q[(at, l2)] * g2
    # ES2 = es2on (onsite, SRULE) + offsite gamma2 (SI Eq 101: 1/[R + 0.5(1/U2A+1/U2B)
    #   exp(-k2x R)], U2 = T32*(1+Gamma*CN), k2x ~ 0.6) + ES3. VALIDATED on F2 (ES3~=0): es2on +
    #   the exact offsite = printed +0.002576 EXACTLY. For HF/H2O the remainder (+0.019/+0.062)
    #   is ES3. The line below is the interim v2 (folded-KO offsite conflated with ES3); the
    #   clean close = es2on + Eq-101 offsite + a SEPARATED ES3. See es2_gamma_blind_diagnosed.
    ours["ES2+3"] = es2on + GE.es_charge_energy(q, zs, B["Rab"], E)
    # ES2 gamma2 DECODED (not fitted -- es2_gamma_gate.py gates it bit-exact vs setespot_'s
    #   own matrix, 8e-11): gamma2 = 1/[R + 0.5(1/U2A+1/U2B) exp(-k2x R)] (Eq 101), with
    #   U2 = T32*ipse*(1+Gamma*CN) [= the MFX first-loop U] and k2x = g2[1] = 0.3300126723.
    #   The full E2 = 0.5 sum q_lA q_lB gamma2 then gates HF +0.17, F2 +0.09, H2O +3.2 mEh vs
    #   printed ES2+3 (so ES3 is small). Earlier T32-only U2 and scanned k2x=0.6 were WRONG;
    #   the dropped ipse factor was the whole error. Productionise: wire E2 with the g-xTB CN
    #   (replacing es2on + the folded-KO), add the small ES3. See es2_gamma_blind_diagnosed.
    # ES3 (SI Eq 129): E3 = (1/6) sum q_lA q_lB (q_A dgamma2/dq_A + q_B dgamma2/dq_B) -- the
    #   CHARGE-derivative of the decoded gamma2. dgamma2/dq = dgamma2/dU2 * dU2/dq; the dU2/dq
    #   kernel is NOT yet pinned (a physical guess dU2/dq=T32*ipse*Gamma gave the WRONG sign and
    #   ~5x too small -- NOT fitted). Needs set3espot_'s charge-intermediates decoded. ES3 is
    #   small (HF/F2 within grade on E2 alone; H2O ~-3 mEh). See es3_decoded_form.
    # ES1 onsite (SI Eq 83b, VERIFIED vs set1espot_ + the SI, validated to sub-mEh):
    #   E1,on = sum_l mu_l * (1 + CN_A*kcn[Z]) * f1(q_A) * q_l
    #   f1(q) = 1 + 0.012*(erf(q-2/3) + erf(q+2/3))    [charge switch; k_s=2/3, k_dis=0.012]
    # mu_l = paramfile shell[5]; kcn = paramfile row0[8]; the erf takes the ATOMIC CHARGE q_A
    # (binary read: HF +-0.412), the linear factor takes CN -- the two were SWAPPED before.
    # Wiring this with the g-xTB covalent CN gates HF -0.32 / H2O -0.23 mEh (was +27.9/+30.5).
    # Interim 0.0165*q_at stand-in below; to productionise: reproduce the g-xTB CN (NOT
    # repulsion.cn_eq47 -- HF 0.306 vs its 0.735) + decode the offsite S-contraction (task #36).
    es1 = 0.0
    for (at, l), qa in q.items():
        z = zs[at]
        es1 += E[z]["MU"][l] * qa * (1 + 0.0165 * qat[at])
    for (at, l), qa in q.items():
        for bt in range(len(zs)):
            if bt == at:
                continue
            zb = zs[bt]
            for lb in range(E[zb]["nsh"]):
                dr = AUFBAU[zb][lb] - E[zb]["ref"][lb]
                gko = 1.0 / (float(B["Rab"][at, bt])
                             + 0.5 * (1.0 / E[zs[at]]["U"][l] + 1.0 / E[zb]["U"][lb]))
                es1 -= dr * gko * qa
    ours["ES1 (charge SIE)"] = es1
    # EHT+ACP electronic piece (no Ex): Tr((H0+A) P)
    tr_part = float(np.sum((B["H0"] + B["A"]) * P))
    # Ex: only H2's gam construction is validated
    if name == "H2":
        gon = 2 * E[1]["cx"] * E[1]["L5"][0]
        gam = np.full((n, n), gon)
        ours["Ex (Mulliken)"] = scf_h2.ex_energy(P, S, gam)
        ours["electronic"] = tr_part + ours["Ex (Mulliken)"]
    ours["atomic core increments"] = sum(INC[z] for z in zs)
    print(f"\n{name}: printed-vs-ours (mEh; + means ours is higher)")
    order = ["electronic", "Ex (Mulliken)", "ES1 (charge SIE)", "ES2+3",
             "ES multipole", "ES total", "Espinpol", "atomic core increments",
             "dispersion", "nuclear repulsion"]
    for t in order:
        if t not in T:
            continue
        if t in ours:
            d = (ours[t] - T[t]) * 1000
            print(f"  {t:24s} printed {T[t]:+12.6f}  ours {ours[t]:+12.6f}  "
                  f"d {d:+9.3f} mEh")
            defects.setdefault(t, []).append(abs(d))
        else:
            print(f"  {t:24s} printed {T[t]:+12.6f}  ours     MISSING")
            defects.setdefault(t + " [MISSING]", []).append(abs(T[t]) * 1000)

print("\n==== THE RANKED DEFECT LIST (worst |d| across systems, mEh) ====")
for t, ds in sorted(defects.items(), key=lambda kv: -max(kv[1])):
    print(f"  {t:36s} worst {max(ds):9.3f}  (n={len(ds)})")
