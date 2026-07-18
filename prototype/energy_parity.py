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
import aes as AES  # noqa: E402
import constants as K  # noqa: E402
import dispersion as DISP  # noqa: E402
import es2_energy  # noqa: E402
import gxtb_engine as GE  # noqa: E402
import mfx  # noqa: E402
import multipole as MP  # noqa: E402
import oracle  # noqa: E402
import overlap as OV  # noqa: E402
import params  # noqa: E402
import repulsion  # noqa: E402
import restart  # noqa: E402

BOHR = K.BOHR
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F"}
P_ = params.parse()


def eeq_charges(raw, n):
    """Parse the binary's printed EEQ(BC) charges (the 'q' column) -- what the repulsion
    Zeff needs (NOT the Mulliken charge; Mulliken gave +3.9 mEh on HF)."""
    lines = raw.split("\n")
    for i, ln in enumerate(lines):
        if "q_CN" in ln and "CN(basis)" in ln:
            out = []
            for ln2 in lines[i + 1:i + 1 + n]:
                m = re.match(r"\s*\d+\s+[A-Za-z]+"
                             r"\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+(-?\d+\.\d+)", ln2)
                if m:
                    out.append(float(m.group(1)))
            return out if len(out) == n else None
    return None

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
    # the CONVERGED state straight from the restart file (P + the printed terms).
    # fock_recon (the old source) only ever supplied its ["state"] -- this same object --
    # but its Fock-inversion assert blocked CH4/NH3 over a matrix the ledger never used.
    st = restart.converged_state(atoms)
    raw = st["raw"]
    T = {k.strip(): float(v) for k, v in
         re.findall(r"^([A-Za-z0-9+() .]+?)\s*:\s*(-?\d+\.\d+)\s*$", raw, re.M)}
    P = st["P"]
    X = np.array(xyz, float)
    B = GE.build(zs, X)
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
    # nuclear repulsion -- Zeff needs the EEQ(BC) charge, NOT Mulliken qat (Mulliken -> +3.9
    # mEh on HF; EEQ -> +0.001). Parse the binary's printed EEQ q; fall back to qat if absent.
    eeq_q = eeq_charges(raw, len(zs)) or [qat[a] for a in range(len(zs))]
    ours["nuclear repulsion"] = repulsion.energy(
        list(zs), X, eeq_q, P_,
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
    # WIRED (was: es2on + GE.es_charge_energy interim, -12.5 mEh): the DECODED gamma2 E2
    # (SI Eq 100-102) with the sqrt/L2 internal CN. Bit-exact vs the binary's own gamma2;
    # computed-CN residual ~0.05 mEh. ES3 (small, +0.024 mEh HF) still to add (task #40).
    ours["ES2+3"] = es2_energy.energy(zs, xyz, q)
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
    # WIRED (was interim 1+0.0165*qat, +27.9 mEh): the full Eq-83b/84 onsite --
    #   mu_l*(1 + k1CN[Z]*CN_A)*f1(q_A)*q_l, CN_A the internal L2/sqrt CN, k1CN=row0[8],
    #   f1(q)=1+k_dis(erf(q-k_s)+erf(q+k_s)) [k_dis=0.012, k_s=2/3], erf takes ATOMIC charge.
    CN1 = es2_energy.coordination(zs, xyz)
    KDIS, KS = 0.012, 2.0 / 3.0
    es1 = 0.0
    for (at, l), qa in q.items():
        z = zs[at]
        kcn1 = es2_energy._P["elements"][z]["row0"][8]
        f1 = 1.0 + KDIS * (math.erf(qat[at] - KS) + math.erf(qat[at] + KS))
        es1 += E[z]["MU"][l] * (1.0 + kcn1 * CN1[at]) * f1 * qa
    for (at, l), qa in q.items():
        for bt in range(len(zs)):
            if bt == at:
                continue
            zb = zs[bt]
            for lb in range(E[zb]["nsh"]):
                dr = AUFBAU[zb][lb] - E[zb]["ref"][lb]
                R = float(B["Rab"][at, bt])   # offsite uses the FULL Eq-101 gamma2 kernel
                gko = 1.0 / (R + 0.5 * (1.0 / es2_energy._U(zs[at], l, CN1[at])
                                        + 1.0 / es2_energy._U(zb, lb, CN1[bt]))
                             * math.exp(-es2_energy.K2X * R))
                es1 -= dr * gko * qa
    ours["ES1 (charge SIE)"] = es1
    # Ex: the DECODED MFX exchange (mfx.py, pass 88) -- setgab_lrao_'s gamma + the 4-index
    # Mulliken energy, x2 closed-shell spin sum; closed to ~1e-8 Eh on all six H..F systems,
    # OFX empirically zero. (Was: the H2-only scf_h2 path, +22 mEh on H2, MISSING elsewhere.)
    gam = mfx.gamma_matrix(zs, X, meta)
    ours["Ex (Mulliken)"] = 2.0 * mfx.ex_energy(P, S, gam)
    # electronic: Tr((H0+A) P) + Ex. The H0 forward assembly is H2-grade only, so the
    # ledger still measures this on H2 alone -- the d there IS the H0 gap, now that Ex
    # no longer contaminates it.
    if name == "H2":
        tr_part = float(np.sum((B["H0"] + B["A"]) * P))
        ours["electronic"] = tr_part + ours["Ex (Mulliken)"]
    # ES multipole: the PORTED AES (aes.py -- SI Eq 116, gdb-extracted erf kernels).
    # CAMM moments from the same restart P; moment integrals in oracle AO order,
    # exactly as aes.run gates it.
    shells_o, _ = OV.build_shells(zs, X, charge=0)
    S_o, D_o, Q_o = MP.moment_matrices(zs, X, shells=shells_o, ao_order="oracle")
    ao_at = np.array([sh["at"] for sh in shells_o for _ in range(2 * sh["l"] + 1)])
    qat_c, dpat, Theta = AES.camm(zs, X, P, S_o, D_o, Q_o, ao_at)
    ours["ES multipole"] = AES.energy(zs, X, qat_c, dpat, Theta)
    # ES total is exactly the sum of the three ES rows (verified against the printout:
    # HF 0.047231+0.045047+0.000011 = 0.092289 = printed) -- a derived row, so its d is
    # the accumulated ES defect (dominated by the missing ES3).
    ours["ES total"] = (ours["ES1 (charge SIE)"] + ours["ES2+3"]
                        + ours["ES multipole"])
    # dispersion: STOCK dftd4 with the derived g-xTB damping (a1=g1[9], s8=g2[9]). This is
    # the revD4 LOWER-BOUND -- it reproduces the binary to ~0.15-0.30 mEh (under grade), the
    # residual being revD4's sigmoidal Mulliken-charge zeta (SI Eq 166-169), not yet ported.
    # So the ledger's dispersion d is the ZETA gap, not an implementation error. See
    # dispersion.py / d4_probe.py.
    d4 = DISP.energy(zs, X)
    if d4 is not None:
        ours["dispersion"] = d4
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
