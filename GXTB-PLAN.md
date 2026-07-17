# g-xTB in this fork: the plan

This is `richapma78/tblite`, branch `gxtb` — **our own implementation of g-xTB inside tblite**,
the same host the method's authors chose (their implementation lives in a *private* tblite fork;
this is the public-side twin). The point of doing it here rather than anywhere else: the day the
authors publish their tblite g-xTB, reconciling is a `git merge`/file diff inside one codebase,
not a rewrite.

## What we have that makes this fully specified

1. **Equations, complete** — the g-xTB preprint SI (ChemRxiv 2025-bjxvt, in the consuming repo's
   `reference/`) derives every energy term AND the analytic gradients.
2. **Parameters, public at v1** — `param/gxtb/` in this repo carries the parameter files from the
   `grimme-lab/g-xtb` **v1.1.0 release** (2025-08-12, the latest public set; sha256 pinned in
   `SHA256SUMS`): `gxtb_parameters` (global + element-wise), `basisq` (the q-vSZP adaptive basis —
   per shell: exponents with a static and a charge-scaling contraction column), `eeq` (the charge
   model). Same license family as this codebase (LGPL-3.0).
3. **A bit-exact oracle** — the standalone v1.1 `gxtb` binary *reads these exact files* (installed
   at `/opt/gxtb-v1` in WSL; verified working: water total −76.43743981 Eh) and prints the energy
   **decomposed by term** (electronic / atomic core increments / dispersion / nuclear repulsion),
   so every term we implement is checkable in isolation, not just the total.
4. **An acceptance instrument** — the ChemRoutes Phase-1 harness (42 frozen-geometry probes, each
   isolating one physical term, with GFN2/g-xTB-v2/DFT labels in a results DB and an
   ordering-gate comparator).

## The terms to implement (SI structure, target file layout)

Mirroring the sibling files here (`src/tblite/xtb/gfn2.f90` etc.) and the layout the private fork's
binary reveals in its strings (`src/tblite/basis/q-vszp.f90`):

- `src/tblite/basis/q-vszp.f90` — the charge-adaptive basis: EEQ(BC) charges computed once per
  geometry, contraction coefficients c = c0 + c1·q, then a standard integral pass.
- `src/tblite/xtb/gxtb.f90` — the Hamiltonian: extended Hückel-type core with the SI's
  distance/CN-dependent scaling, **fourth-order** charge self-consistency (vs GFN2's third),
  Mulliken-approximated range-separated exchange, atomic correction potentials (ACPs; f-projectors
  on the heavy blocks), spin polarization constants, atomic core increments (these put totals on
  the all-electron wB97M-V/def2-TZVPPD absolute scale — the v1 oracle prints them separately).
- dispersion: revD4 (upstream tblite already carries D4 machinery to extend).

## Gates (declared before the build, per the consuming repo's discipline)

- **G1**: parse `param/gxtb/*` and reproduce the v1.1 oracle's *term decomposition* on water to
  1e-6 Eh per term.
- **G2**: reproduce the oracle's totals across the 42 harness probes (bit-tight target; declare
  the tolerance before running, fail = find the bug, never tune).
- **G3**: only after G2 — distill the parameter file toward the **v2** binary on harness probes
  (start from v1 values; this is fitting, not reverse engineering). The harness ordering gates
  (acyl ladder, Pd ladder) must stay MATCH throughout.
- **Merge day**: when upstream publishes g-xTB, diff/merge their implementation against this
  branch; whichever survives, the harness gates decide equivalence.

## Consumers

- The ChemRoutes sandbox dispatcher (`sandbox/qm.py`) already has the seam: its tblite backend
  takes the method as a string — this fork's build adding `"g-xTB"` makes it available with no
  caller changes.
- xtb consumes tblite as a meson subproject (`subprojects/tblite.wrap`); pointing the wrap at this
  fork yields a full xtb binary carrying our g-xTB, exactly how the authors' distribution works.
- A later GPU port (PySCF/CuPy batch-across-candidates) validates against this implementation via
  the same harness.

## Status log

### 2026-07-16 — G1 groundwork: the oracle is decoded end-to-end, three instruments live

- **`prototype/` is the build's first stage** (decided): the physics is proven in transparent
  Python against the oracle FIRST; the Fortran port into `src/tblite/` happens once G2 passes.
  Porting proven physics is mechanical; debugging physics in Fortran is not.
- **`prototype/oracle.py`** — drives the v1.1 binary and parses EVERY printed intermediate:
  EEQ(BC) charges/CNs, the q-vSZP per-atom adaptation, basis dims, SCF trace, eigenvalues, shell
  populations, WBO, Janak IP/EA, and the full TERM DECOMPOSITION. `-grad` gives numerical forces
  (the reference for G4). Selftest: water, 10 terms parsed.
- **Term 1 DONE — atomic core increments**: pure per-element constants, read off single-atom
  runs (`prototype/data/increments-v1.json`, 13 elements). incr(H) = 0 (no core). **Additivity
  VERIFIED to 1e-7** on water/HCl/CO. The values mirror def2-ECP conventions (Br all-electron
  −2569.3; I effective-core −293.6; Pd −123.9) — as expected for a method targeting def2-TZVPPD.
- **`prototype/basisq.py`** — the q-vSZP file layout is decoded and pinned: per element, a header
  (Z + the three adaptation coefficients of `q + a·q² + b·CN^0.5 + c·q·CN`) then shells of
  primitives with TWO contraction columns (static c0, charge-scaling c1). All 103 elements parse.
  Conventions measured on lone atoms: the oracle's `npr` counts primitives × Cartesian components
  (O: 24 = 6s + 6p×3); `acpsao` is 16 PER ATOM = 1+3+5+7 — **the ACPs carry s,p,d,f projectors on
  every element**.
- **`prototype/eeq.py`** — EEQ(BC) file layout decoded: a date-stamp header (same vintage stamp
  the binary prints) + 103 rows × 10 floats (8 meaningful, 2 zero-padded). No globals in the file.
- **42 harness-probe labels banked** from the v1 oracle (ChemRoutes `sandbox/harness.py --teach
  gxtb-v1`, 0 failures incl. charged/open-shell probes via `.CHRG`/`.UHF`) — gate G2's target set
  now exists as data.

### Gate G4 (added): analytic gradients

The v1 binary is numerical-only; the SI derives the ANALYTIC gradients completely; the v2 binary
ships them. Our build implements the SI's analytic derivatives after G2, gated three ways: match
our own finite differences (internal consistency), match the v1 oracle's `-grad` numerical forces
on the probe set (external truth at v1 parameters), and only then time against the v2 binary.
This also restores what the ecosystem note observed ("initial releases lack analytical gradients,
making optimizations slow") — with gradients, optimization/frequencies at g-xTB cost stop being
the bottleneck.

### Solvation strategy (the layer upstream will not ship soon)

- **Borrowed-cavity ansatz — GATE FAILED, recorded, not tuned** (ChemRoutes
  `sandbox/validate_solvation.py`): E_gxtb(gas) + [GFN2-ALPB shift] reproduced the wB97M-V+SMD
  acetonitrile ladder's separated structure (acetate ≫ chloride > heavy cluster) but inverted two
  rungs INSIDE the referee's 0.2 kcal/mol degeneracy. Strict full-ordering gate: FAIL → solvated
  questions stay with the DFT referee. Post-mortem: the neutral-neutral acyl set is a WEAK
  solvation test (rung effects −1.9..+0.5 kcal/mol).
- **The native plan**: parameterize ALPB (tblite already carries the solvation machinery) for
  g-xTB's charges, trained against the GPU-DFT+SMD labels the ChemRoutes escalation tier banks as
  a side effect of normal work. Gate to declare before fitting: ordering + magnitudes on a
  **charged-species** set (acylpyridinium-type adducts — where solvation is tens of kcal/mol and
  the engine's catalysis screen actually lives), not the neutral set that just failed weakly.

### Codebase audit verdicts (2026-07-16, full read of this fork — details in ChemRoutes docs/engine-map.md)

**Eases (much of the scaffolding already exists upstream):** h0spec hooks `get_q1shift/get_q2shift`
(+`kq1/kq2` storage, `dsedq`), diatomic-frame σ/π/δ scaling end-to-end (`get_diat_scale`,
`integral/diat_trafo.f90`, `overlap_diat`, working consumer in CEH), EN-weighted CN, complete spin
polarization (needs only g-xTB W constants), `multicharge::get_eeqbc_charges` (the exact EEQ(BC)
q-vSZP needs, already a guess), CEH proves f-shells + Z≤103 through basis/integrals,
self-consistent D4 via the dftd4 subproject, and a generic `interactions` container for bolt-on
terms. 4th-order charges = copy `coulomb/thirdorder.f90`. Atomic increments = trivial classical
container (data in prototype/data/).

**Hard parts (the real work, in order):** (1) q-vSZP — per-ATOM charge-adaptive contractions
rebuilt each geometry (+ ∂c/∂q), no existing path; (2) charge-dependent H diagonal into the SCC
loop (self-energy currently computed once pre-SCF); (3) Mulliken range-separated exchange — needs
a new off-diagonal Fock pathway, nothing similar exists; (4) ACP projectors — absent; (5) custom
parameter parser (the TOML `element_record` lacks every g-xTB field — mirror gfn2.f90's
module-owned parameters instead).

**Method registration** touches ~5 dispatch sites: `app/driver_run.f90`, `app/driver_param.f90`,
`src/tblite/api/calculator.f90` (+ header), `python/tblite/interface.py::_loader`, plus
meson/CMake lists.

**Solvation layer seam confirmed LOW-RISK:** ALPB/GBSA/CDS parameters are selected BY METHOD
STRING (`solvation/data/*` + `get_alpb_param(..., method, ...)`) and consume `wfn%qat`, which
g-xTB produces — the work is new `"gxtb"` parameter tables (fit against ChemRoutes' banked
DFT+SMD labels) and the pre-declared charged-species gate, not plumbing.

**No numerical-gradient fallback exists in the library core** — every term ships analytic
gradients (FD lives only in unit tests as verification). The SI's gradient derivations are
load-bearing for G4.

### 2026-07-16 (second push) — Term 2 GATED: the EEQ(BC) charges are OURS now

- **Identification, exhaustive**: the v1 `eeq` file is the PUBLISHED eeqbc2025 parameter set
  (Froitzheim/Mueller/Hansen/Grimme, J. Chem. Phys. 2025, 162, 214109) — byte-identical across
  all 103 elements × 8 columns (`prototype/extract_tables.py` re-proves it on every run; column
  order chi, eta, rad, kcnchi, cov_radii(raw; working radius = half), kqeta, kqchi, cap). The
  reference implementation is multicharge's `model/eeqbc.f90` + factory `new_eeqbc2025_model`
  (constants: kcn 2.0, norm_exp 0.75, kbc 0.60, kcnrad 0.14, cutoff 25 Bohr; EN = Pauling/3.98
  with actinide patches; avg_cn from the published array — not in the v1 file).
- **Implementation**: `prototype/eeqbc.py` — erf-CN, EN-weighted local charge, bond-capacitance
  Maxwell matrix, CN-scaled Gaussian widths, bordered constrained solve. One faithful quirk:
  the pair vdW radii enter as their ÅNGSTRÖM literals against distances in BOHR (multicharge
  converts to au at declaration and back at use) — reproduced as the code runs, and the gate
  confirms the oracle's binary does the same.
- **GATE PASSED on the full probe set**: water, HCl, CO, NH4+ (charged), AcCl, PdCl2 (transition
  metal) — CN, q_loc (the oracle's q_CN column) and final q all within the oracle's 4-decimal
  print precision (max |diff| ≈ 5e-5).
- New sibling reference clones: `C:\Projects\multicharge` (eeqbc reference + published params)
  and `C:\Projects\mctc-lib` (CN counting functions + vdW/EN data tables), both read into
  `prototype/data/mctc-tables.json`.
- **Next in the chain**: CN(basis) — the q-vSZP CN parameterization (the SI's "four different CN
  parametrizations"; the oracle prints CN(basis) per atom = direct gate) → q_eff assembly
  (`q + a·q² + b·CN^0.5 + c·q·CN`; the oracle's AO-setup block prints every piece per atom) →
  overlap integrals in the adapted basis.

### 2026-07-16 (third push) — Term 3 GATED: the basis CN and q_eff; the adaptive basis is ours

- **The basis CN was hiding in plain sight**: not in the SI (Eq. 47 is the HAMILTONIAN CN, a
  different one of the four) and not CEH's CN (tested and refuted: CEH radii give pair count
  0.86 where the oracle says 0.41). The authoritative source is the q-vSZP SETUP TOOL
  (`grimme-lab/qvSZP`, new sibling clone) — `ncoord_basq`: erf count with **kn = −3.75** over
  rc = SUM of Pyykkö–Atsumi 2009 covalent radii (metals −10%), Å→Bohr. A first pure-fit attempt
  from oracle diatomic scans (`prototype/fit_basis_cn.py`, kept as the record of the method)
  bracketed the constants; the tool's source pinned them.
- **q_eff was decoded directly against the oracle's AO block before reading anything**: with the
  basisq element-header triple (h1,h2,h3), q_eff = (q − h2·q²) + h1·√CN + h3·q·CN reproduced
  every printed component (dq/dcn/dqcn/total) to all decimals — then SI Eq. 28 confirmed the
  form (its `a` is −h2; k0 = 1 default). `prototype/adapt.py` implements both.
- **GATE PASSED on all six probes** (water, HCl, CO, NH4+, AcCl, PdCl2): CN(basis) ≤ 4e-5,
  q_eff ≤ 5e-6 vs oracle.
- Consequence: c = c0 + c1·q_eff per primitive is now fully computable by our code for ANY
  geometry — the charge-adaptive basis (the audit's "single biggest structural lift" for the
  Fortran port) has a working, gated reference implementation.
- Note for the Fortran port: SI Sec. 1.2 also defines a SECOND, scaled basis variant used only
  for the EHT Hamiltonian (scaled exponents + scaled k0/k2/k3, element-wise, from the main
  parameter file) — "decouples the overlap from the effective Hamiltonian". Two basis builds
  per calculation, one q_eff.
- **Next in the chain**: overlap integrals in the adapted basis (oracle gates: nsao dims ✓
  already, then eigenvalue/population checks once H0 exists); the diatomic-frame scaled overlap
  (SI Sec. 1.3) machinery already exists in this fork (`integral/diat_trafo.f90`).

### 2026-07-16 (fourth push) — Term 4 GATED: the overlap engine, machine-exact vs PySCF

- **`prototype/overlap.py`**: Obara–Saika cartesian primitive overlaps (general l), a
  cartesian→real-spherical transform CONSTRUCTED NUMERICALLY (least-squares against scipy's
  spherical harmonics on random unit vectors — exact for pure-l polynomials, no transcribed
  tables to mistype; handles both scipy generations' sph_harm APIs), PySCF's normalization
  convention mirrored (axis-normalized primitives + unit-norm contracted shells).
- **The one real bug found by the gate**: PySCF orders p shells (x, y, z), not m = −1..+1 —
  a permutation worth 0.64 in max|dS| until diagnosed against `gto.mole.cart2sph` directly
  (which also revealed pyscf's per-l scalar, irrelevant post-normalization). One row reorder
  → machine precision.
- **GATE PASSED**: max|dS| ≤ 5.5e-15 on water (s/p), AcCl (d), PdCl2 (d on a TM), and CeO
  (**f(7) on cerium — lmax 3 covered**), all in the ADAPTED basis built by the gated chain
  (eeqbc charges → basis CN → q_eff → c = c0 + c1·q_eff).
- Note: the ORACLE's own AO normalization convention is deliberately not asserted by this gate
  (nothing printed exposes S directly); it gets pinned at the H0 stage where eigenvalues and
  Mulliken shell populations become observables. The gxtbrestart file (binary MOs?) remains an
  unexplored shortcut if eigenvalue-stage debugging ever needs the density directly.
- **Next in the chain**: the SECOND basis variant (Hamiltonian-scaled exponents + scaled
  k0/k2/k3 from the main parameter file — SI Sec. 1.2 end), the diatomic-frame scaled overlap
  (SI 1.3; machinery exists in this fork), then H0 (SI 1.7) with eigenvalue gates. That stage
  requires decoding the main `gxtb_parameters` file layout — the largest remaining decode.

### 2026-07-16 (fifth push) — Term 5 GATED (repulsion); the parameter file is CRACKED OPEN

- **The decisive new instrument: perturbation mapping** (`prototype/perturb_map.py`). The oracle
  re-reads `~/.gxtb` every run and prints every term — so nudge ONE float, rerun fixed probes,
  record which observables moved. One sweep named ~40 slots mechanically
  (`data/perturb_map.json`): repulsion element params (α0, R0, Zeff0, kCN, kq | L1[0..4]), the
  shared Eq.47 CN radius (L1[5] — moves electronic AND repulsion together), penetration globals
  (G1[3]=kpen1, G1[8]=kpen1_H/He, G2[2,3,4]=kpen2..4 — sensitivity decaying with power,
  textbook), the Eq.47 steepness (G2[0], stored NEGATIVE — resolving the SI's apparent sign
  typo), dispersion globals (G1[9], G2[9]), multipole globals (G2[6,7]), per-shell rows
  L2..L7 = levels / kCN-level / ζ-scale-ish / exchange / Hubbard γ / first-order, L9 = ACPs.
  Two long-range/charged discriminator probes split α0 from kCN (kCN goes dark at CN→0) and
  located kq2 at L8[6]. **L8[5] moved the increments line by exactly N·δ — Term 1's parameter
  found in situ**, and the file's 79 increments match Term 1's oracle measurements at 4.5e-9
  (`params.py` selftest re-proves it).
- **File map**: 20 globals + 79 element blocks (H..U; 4f/5f gaps: 59–69, 90–91, 93–103 absent),
  each block L1(10 scalars) + L2..L7(4-wide per-shell) + L8,L9(8-wide). `prototype/params.py`
  parses it with names where pinned.
- **Term 5 — repulsion (SI 1.6) implemented and GATED** (`prototype/repulsion.py`): worst
  |dE| = 4.75e-9 Eh over 13 cases (full H2 curve, HCl, NH4+ charged, water, AcCl, PdCl2).
  Three SI-extraction ambiguities settled EMPIRICALLY, each recorded in the code: (1) the
  exponent offset is (R + R0)^1.5; (2) rc in Eq.47 is the arithmetic MEAN (as the SI says in
  words); (3) α_AB is the HARMONIC mean 2αAαB/(αA+αB) — Eq. 55's extraction lost a factor 2
  (the H2 tail refuted the halved form by 150×); (4) BONUS: Zeff = zeff0·(1 − kq·q − kq2·q²) —
  the quadratic sign was measured by FD probing (oracle/model ratio −1.000 exactly), and the
  linear term validated at ratio 1.000 with q = the EEQ(BC) charge.
- Chain now: increments ✓ charges ✓ adaptive basis ✓ overlap ✓ repulsion ✓. Classical energy
  complete except dispersion (revD4 — a parameter selection on the dftd4 library, its two
  globals already located). **Next: the EHT Hamiltonian** — levels/CN shifts (L2, L3 named by
  the map), the Ham-basis scalings, diatomic-frame overlap, eigenvalue gates.

### 2026-07-16 (sixth push) — the restart file cracked: the converged state is now an INPUT

- **Why it matters**: with the oracle's converged density matrix P in hand, every electronic
  term (ES1, ES2+3, multipole ES, Mulliken exchange, spin) can be gated INDEPENDENTLY against
  the printed per-term decomposition — before our own SCF loop exists. The SCF then becomes
  fixed-point iteration around already-gated pieces. `prototype/restart.py`.
- **Format** (measured): Fortran sequential records — record 1 = packed-triangular symmetric
  DENSITY (row-lower packing), record 2 = the MO matrix (column-major; big systems store only
  occupied columns — the matrix goes singular, expected).
- **The S-extraction trick**: MO orthonormality CᵀSC = 1 means S_oracle = C⁻ᵀC⁻¹ — the oracle's
  TRUE overlap matrix recovered from the restart alone. Water: S_oracle equals OUR overlap to
  1.9e-8 — the oracle's s/p basis, normalization and ordering are exactly ours.
- **The d-ordering measured**: HCl showed a pure two-entry swap; a tilted-HCl probe (all five
  d-components lit) solved the full signed permutation at residual 1.1e-8, signs all +1:
  **oracle d order = (x²−y², z², xy, xz, yz)** vs PySCF's m = −2..+2. Wired into
  `overlap.py` as ao_order="oracle" (f ordering unmeasured — refuses, by design, until a
  lanthanide probe measures it).
- **RESTART GATE PASSED**: parsed P + our oracle-ordered S reproduce the printed Mulliken shell
  populations on water/HCl/AcCl (worst 4.9e-4 = print precision; Tr(PS) = nel to 1e-4). The
  PySCF machine-precision overlap gate is untouched and still passes.
- **Next**: with P as input, implement + gate the electronic terms one at a time against the
  printed decomposition — ES2+3 (isotropic 2nd/3rd order; Hubbard γ slots already named at L6),
  ES1 (first-order, L7), multipole ES (globals G2[6,7]), Mulliken exchange (L5 + G1 slots),
  spin (Espinpol); then H0 closes the electronic energy via Tr(PH0) and the SCF loop follows.

### 2026-07-16 (seventh push) — first-order DECODED (fractional references!); onsite 2nd-order kernel CORNERED

- **The fractional-reference discovery**: lone neutral atoms converge to aufbau shell populations
  yet print NONZERO ES1/ES2+3. Perturbing oxygen's L7 slots moved ES1 by ±0.325 — equal and
  opposite for s and p — exposing FRACTIONAL reference occupations. SI Tab. 1 then confirmed:
  the references are averaged wB97M-V/q-vSZP Mulliken occupations (O: s 1.673 p 4.326; the FD
  measured 1.675/4.325). Measured for 8 elements incl. Si's three shells — all match the table
  to ~4e-4, pinning the convention **q_l = refocc − pop** (the sign opposite Eq. 99 as printed).
  ES1 is exactly LINEAR in the L7 chemical potentials, so the FD extraction is exact — refoccs
  for elements beyond the table's Z≤18 are measurable the same way (needed: the file does NOT
  store them; no element slot other than L7/L6/L5/L1[9] moves the atom's terms).
- **ES1 (atoms) = Σ_l µ⁰_l·q_l with µ⁰ = the L7 row** — structure verified (switching f(0)=1).
  Molecule pieces still to assign: k^(1),CN (element), k_dis/k_x/k_s (globals), offsite via γ⁽²⁾
  with Δρ⁰ = aufbau − refocc.
- **Onsite second-order: the kernel is NOT any standard mean.** With q measured and U = L6, all
  of harmonic/geometric/arithmetic fail with element-VARYING ratios (0.59–2.5). Three exact
  constraints for O (energy + two FD slopes, globals PROVEN inert for atoms — a full 20-slot
  sweep moved nothing) give γ_sp(1.4936, 0.8473) = 0.5321, ∂γ/∂U_s = 0.3608, ∂γ/∂U_p = 0.6814;
  Euler's test on these says γ_sp is HOMOGENEOUS OF DEGREE ≈ 2.1 in the Hubbards — not the
  SI's Eq. 101 R→0 limit. Next experiment (designed, not yet run): SET the L6 values to chosen
  grids (equal-U scan for the diagonal identity, magnitude scan for homogeneity degree,
  asymmetry scan for the mean's form) and chart γ(U_s,U_p) directly — the identification is
  2-D function mapping, the instrument exists (`es_atoms.py` carries the measurement machinery;
  its declared gate currently FAILS and stays failing until the kernel is measured, not tuned).

### 2026-07-16 (eighth push) — the onsite 2nd-order kernel MEASURED; the private fork's file tree read off the library

- **The user's hint paid off — the distributions carry structure**: `ar t libxtb.a` on the Linux
  install lists ALL 710 object files of the build INCLUDING the private tblite fork's source
  tree by name: `basis_q-vszp.f90` (as predicted), `classical_increment.f90`,
  `coulomb_multipole_gxtb.f90` (a g-xTB-specific multipole), and the decisive
  **`coulomb_thirdorder_onsite.f90` + `_twobody.f90`** — onsite third-order exists as its own
  container, exactly what the atom anomalies demanded. Also: the "Windows" release asset is a
  REAL native g-xTB build (`--gxtb` works; inner folder confusingly named xtb-bleed-windows),
  extracted at ChemRoutes reference/gxtb-win with gcc .mod interface files + libxtb.a.
- **The survey experiment worked exactly as designed** (SET the L6 Hubbards, don't nudge):
  equal-U scan → ES2+3(O atom) is DEAD CONSTANT (C0 = 0.00012001, U-independent) — the E2 part
  vanishes with the near-zero total charge and the leftover is a pure higher-order onsite
  constant. Asymmetry scan (U_p = 1, U_s ∈ 0.5..3) → the U-dependent part is EXACTLY
  K·(U_s−U_p)²/(U_s+U_p) (four scan points, ratios constant to 4 digits; even in ΔU to <0.2%,
  so the odd-in-ΔU third-order shape is ABSENT here) — and that is precisely the
  harmonic-mean-kernel form: E2 = s·[½ Σ q_l q_l' γ_harm(U_l,U_l')]. The measured K reproduces
  the earlier FD slopes to 3 digits (dES23/dU_s +0.01493 vs measured +0.0150).
- **The scale s is per-ELEMENT, not global** (measured: C 0.417, N 0.496, O 0.582, S 0.364,
  Cl 0.380) and decreases smoothly with |q| for the 2-shell elements — a charge-damping of the
  second-order onsite, closed form not yet identified (the 20 globals are PROVEN inert for
  atoms; s is not stored per element in any slot that moves ES23 — likely a hardcoded damping
  function of the charges, like avg_cn was hardcoded). C0 measured per element (C 0.00909,
  N 0.00161, O 0.00012, S 2.1e-5, Cl 1.4e-6) — roughly quartic in the shell charges,
  the fourth-order onsite candidate.
- **Next experiment (designed)**: charged atoms (.CHRG −1/+1/+2) sweep q_l systematically at
  fixed U → chart s(q) and C0(q) → identify the damping's closed form and the higher-order
  onsite polynomial; then molecules (offsite kernel with its exp(−k(2),x·R) screening).

### 2026-07-16 (ninth push) — the charge scans: f(q) charted, E3 measured cubic, hardcoded-constants PROVEN

- **The first-order switching f(q) is measured and matches the SI's form** — as an ODD deviation
  (my earlier even-parity reading of Eq. 83b was wrong; the erf-pair SUM is odd): FD products on
  the L7 potentials across O charge states −2..+2 give f = (0.9764, 0.9831, 1, 1.0165, 1.0234),
  saturating erf-pair shape; constants k_dis ≈ 0.012–0.013 with saturation near |q|~2 (exact fit
  deferred to molecular fractional-charge data). The FD also re-measured q_l per ion, confirming
  aufbau populations and fractional references throughout.
- **The onsite third-order is ALIVE and exactly cubic in the ATOMIC charge**: equal-U ion data
  splits into even/odd parts; the odd part is 0.0418·q_A³ (±1 vs ±2 coefficients agree to 1%) —
  at U=1 this pins the Γ·τ combination at ≈ 0.251 for O. The even part at ±1/±2 scales as q_A²
  (ratio 4.21 ≈ 4), i.e. genuine second-order of the now-charged atom with s ≈ 0.60–0.64 —
  consistent with the neutral measurement's s(O) = 0.582 drifting with charge.
- **PROOF: all onsite-electronic constants are hardcoded in the binary.** Two null sweeps at
  O²⁺ — all 20 globals, then all 50 O-block slots — moved ES1/ES23/Ex/spin through NOTHING but
  L5 (exchange), L6 (Hubbards), L7 (chemical potentials), L1[9] (spin W). The f-constants, the
  third-order rule, the s(q) damping, and the fourth-order onsite all live in code, not in the
  parameter file (as avg_cn did). The 20 file globals serve MOLECULAR physics only (repulsion
  penetration, CN steepness, dispersion, multipole, and the yet-unassigned offsite
  kernel/exchange ranges).
- **Consequence for strategy**: the onsite functional forms must be charted (more charge states,
  more elements, fractional charges via molecules) or lifted from procedure signatures — the
  private tblite .mod files do NOT ship (checked; only xtb-side delegation interfaces), so
  charting remains the instrument. The offsite/molecular terms, by contrast, have their
  parameters IN THE FILE and are the natural next implementation targets alongside the
  remaining onsite charting.

### 2026-07-16 (tenth push) — HOUSE RULE: no method constant lives in code

Per Richard: everything decoded or measured goes into a TABLE, loaded once at start, like the
official parameter files. Implemented: `prototype/data/derived-constants.json` (one home for
every constant NOT in gxtb_parameters/basisq/eeq — basis-CN radii+kn, EEQ(BC) factory constants
incl. the unit quirk, repulsion structural readings, the fractional reference occupations,
oracle conventions incl. the d-permutation and restart format, and the measured-but-not-closed
onsite functions as explicit calibration targets with provenance per entry) +
`prototype/constants.py` (the single loader). All five consuming modules refactored; the FULL
gate suite re-run and green (eeqbc, adapt, repulsion 4.75e-9, overlap machine-precision,
restart Mulliken) — the refactor is proven behavior-identical. The Fortran port inherits this:
our tblite implementation will LOAD what upstream hardcodes, making the fork strictly more
configurable than the binary it reproduces.

### 2026-07-16 (eleventh push) — the onsite charting campaign: the s-RULE FOUND

- **The clean channel**: fitting each ion's equal-U scan as A + B·u + C·u² decomposes the onsite
  electronics exactly — B = (s/2)·q_A² (even, quadratic to 5 digits across ±1/±2), C = k3·q_A³·u²
  (cubic to 4 digits) — while A absorbs the ion shell-redistribution mess, which is why earlier
  neutral-only estimates of s were contaminated.
- **s HAS A FORMULA**: s = a_p·Nval + b_p, LINEAR in the valence-electron count within a period
  (9 elements, max deviation 2e-4): period 2 (a 0.08247, b 0.0920 — Li at Nval=1 lands EXACTLY on
  the line through C/N/O/F), period 3 (a 0.04193, b 0.1378 — Si/P/S/Cl), H 0.4726, Br ~0.437
  (period-4 anchor). The slope HALVES from period 2 to 3. A hardcoded closed form, now charted;
  period 4/5 completion + the (a_p, b_p) generating rule remain.
- **The third-order carries a U² rule**: the equal-U u² coefficient = k3·q_A³ exactly (C, O all
  four ion states); k3 measured for 8 elements — row-2 C→N→O linear, F anomalous (needs the
  per-state check; d-block extraction needs population tracking — Ge showed ±1 states disagreeing
  0.28 vs 0.40 through shell shuffling).
- All banked in `data/derived-constants.json` (s_rule, third_order_u2_rule, raw scans in
  `data/onsite-charts.json`); constants loader green.

### 2026-07-16 (twelfth push) — campaign v2: the third-order SHAPE identified; F anomaly dissolved

- Per-point POPULATIONS (v1's missing ingredient) + the s-rule subtracted pointwise isolate the
  third-plus-higher-order residual on a 9-point (U_s,U_p) grid per ion state.
- **The onsite third-order kernel is the DEGREE-2 form**: E3 = c1·q_A·Σ_ll' q_l q_l'(U_l+U_l')²/4
  — the Eq. 132 (U+U')²/4 structure applied onsite (the SI's onsite line reads 'harmonic' but its
  extraction was garbled; the measured u² law and now the full 2-D shape both say degree 2).
  c1 per element: C 0.0197, O 0.0318, F 0.0163, S 0.0075 — consistent across ± ion states
  (F: 4-digit agreement), C⁻ captured at 1e-3 over a 0.17 Eh range.
- **The F 'anomaly' dissolved**: with the right shape, fluorine was correct all along; the
  apparent row-2 linear k3 trend was the artifact.
- Open: C⁺/O²⁺ carry unmodeled (shell-resolved quartic?) structure; c1's element rule unknown;
  raw grids in data/onsite-charts2.json. NOTE the strategic unlock: for NEUTRAL molecules the
  atomic q_A are small (|q|≲0.3 → E3 ~ 1e-3, E4 smaller), so the molecular OFFSITE gate can
  proceed with the charted onsite pieces before the ion-regime forms fully close.

### 2026-07-16 (thirteenth push) — the offsite front opened: kernel Coulombic tail confirmed; the CN-slots channel found

- **HCl stretch (10 points, populations per point)** with the charted onsite pieces subtracted:
  the offsite ES2+3 tail is CLEAN Coulomb (q·q'/R to ~10% at R=8) and the Klopman-Ohno form is
  roughly right unscreened; the short-range mismatch is the U CN-dependence (kU slots), not yet
  fitted. Data in `data/offsite-hcl.json`.
- **A hypothesis executed and killed properly**: the offsite ES1's overlap-like fast decay
  (factor 28 over R=2..3.5 while the kernel drops 1.7x) suggested an overlap-induced reference
  shift -- computed with our gated overlap engine, it is IDENTICALLY zero (Mulliken of an
  atom-diagonal reference density equals refocc regardless of S; the math nulls it). The real
  source: the CN-dependence of the onsite chemical potentials mu = mu0*(1 + k1CN*CN) -- CN(R)
  decays erf-fast, exactly the measured shape, and the atom-derived onsite model had CN = 0.
- **Consequence**: the stretch data IS the measurement channel for the per-element CN slots
  (k1CN, kU) plus the global k2x -- a 5-parameter joint fit against 20 observables, with the
  Eq. 47 CN computable exactly from the repulsion-era decode. Queued next with fresh budget;
  candidate element slots L1[6]/L1[7]/L1[8] to be value-matched then FD-verified.

### 2026-07-16 (fourteenth push) — molecular slots ASSIGNED by response shapes; a bad fit rejected properly

- A 5-parameter forward fit on the stretch converged numerically (rms 6.5e-4) at k1CN ≈ −23 —
  one hundred times the SI's stated scale. REJECTED as wrong-structure (smooth curves + free
  parameters hide missing physics). The honest instrument instead: FD-probe each candidate slot
  at R = 2/3/5 and read its ROLE off the response-vs-R shape.
- **Assignments (shape-identified, file values sane):** L1[8] = k1CN (µ CN-dependence; ES1-heavy,
  erf-fast decay, H/Cl signs opposite; the file's H value 0.775 reproduces the measured stretch
  signal magnitude −0.011 vs −0.019 with Cl's share pending); L1[6] = kU (Hubbard
  CN-dependence); L1[7] = the Hamiltonian-basis scale (the only L1 slot alive at dissociation —
  it moves the CHARGES); G1[0] = short-range exchange global; G1[5] = R-independent exchange
  component (range-separation signature); G2[1] = offsite-kernel global (k2x candidate).
- Next: the VERIFICATION forward model — onsite(charted) + CN terms with FILE values + offsite
  kernel — against the full stretch; the remaining structural unknown is Eq. 86's Δρ (the
  garbled extraction; aufbau−refocc gives a near-cancelling sum that undershoots 30×).

### 2026-07-16 (fifteenth push) — Eq. 85 read VISUALLY: the paper itself has the typo; first order CLOSED

Rendered SI pages 47–48 to images (pdftoppm) and read the typeset math directly: Eq. 85 is
printed with inconsistent subscripts (Δρ⁽¹⁾_{0,l_B} = ρ_{0,l_A} + Z_{l_A}) — a typo IN THE
PAPER, resolved by its prose: Δρ = aufbau − refocc per shell, static and small. Combined with
the previous push's k1CN identification, the first-order term is fully accounted for (µ = L7,
k1CN = L1[8], switching globals hardcoded, offsite genuinely tiny). Method note: when an
extracted equation resists sense-making, render the page and READ it — the text layer and even
the print can lie in different ways.

### 2026-07-16 (sixteenth push) — Eqs. 85/86 CANONICAL (user's screenshots); stretch tail VERIFIED with file values

- Richard supplied zoomed screenshots of Eqs. 85/86 (ChemRoutes reference/SI_Equation85.png,
  SI_Equation86.png) — now the canonical reference. Confirmed: the subscript typo is in the
  paper; Δρ = aufbau − refocc per shell; Eq. 86 exactly as implemented; the parameter census
  (µ shell-resolved, k1CN element-wise, three switching globals) verbatim.
- **The verification forward model (FILE values only)**: HCl stretch, both observables. TAIL
  (R ≥ 4) VERIFIED — ES2+3 within 4–9e-4, ES1 within 3e-4: the offsite kernel, the static Δρ,
  and the asymptotics are all correct as decoded. SHORT RANGE carries a smooth CN-shaped
  residual (~1e-2 at R=2 on both observables) and an FD cross-check quantifies the miss:
  measured dES1/dk1CN(H) = −0.0445/unit vs the model's −0.0179 — factor 2.5 — so the µ-CN
  coupling's STRUCTURE (which CN convention / what it multiplies) needs one more decode round;
  the slot itself (L1[8]) stands. Same suspicion applies to kU (L1[6]; H's value 0.81 is odd
  against the SI's 'typically <0.2'). Next: FD-shape decode of the CN-term structure (rc
  sum-vs-mean, per-shell vs summed, f-coupling) at fixed geometry.

### 2026-07-16 (seventeenth push) — Richard's text extraction: three finds

- `ChemRoutes reference/Supplement_Information_TextCopyOnly.txt` (direct copy-paste, passes the
  no-biology gate, stays in place). Three yields: (1) **one internal CN** for
  repulsion/Hamiltonian/TB1+2/exchange (SI 1.4 verbatim — the convention our repulsion gate
  already verified; Eq. 84's CN is settled); (2) **the exchange kernel's form surfaced**
  (Eq. 149: η + (1−η)·erf(ωR) — textbook range separation, matching G1[5]'s R-independent FD
  signature); (3) **Eq. 55 prints the unhalved combination in BOTH extractions** — the paper and
  the v1.1 binary genuinely differ by the factor 2 our gate measured; the binary is the
  authority for this build.
- The 2.5× µ-CN slope gap now has a prime suspect: DENSITY RELAXATION in the FD response
  (the printed-term derivative includes the SCF charge shift; the direct term doesn't) — next
  round measures dq/dslot from the printed populations under FD to separate direct from
  relaxation, and probes whether the hardcoded onsite s/c1 carry CN-dependence (charted at
  CN=0 only so far).

### 2026-07-16 (eighteenth push) — relaxation RULED OUT; the µ-CN coupling shape is nonstandard

- FD on k1CN(H) at five distances WITH population readout: the density-relaxation contribution
  to dES1/dk1 is ~3% at R=2 (computed from the measured dq under perturbation) — **the 2.5×
  slope gap is NOT relaxation; it is the direct coupling's structure.**
- The extracted direct shape, direct/(µ·q) over R = (0.557, 0.282, 0.097, 0.022), follows
  NEITHER Eq.47-CN convention (mean-rc decays too fast: 0.219→0.0005; sum-rc too slow:
  0.816→0.361) and is not a single erf of any (k, rc) (fits miss by ~25%). Either the µ-CN
  uses its own counting parameters, or the coupling multiplies something beyond µ⁰·q.
- Instrument limitation identified: printed populations carry 4 decimals → dq under FD is
  noise-limited (the R=4 row is unusable). UPGRADE DESIGNED: charges from the restart DENSITY
  via our gated overlap (machine precision) for all FD work — serves every future round.
- Data: this round's slopes recorded here; raw in the transcript. Next: restart-precision FD
  repeat + shape-test against per-element-steepness variants and µ-independent couplings.

### 2026-07-16 (nineteenth push) — the µ-CN coupling is √CN (machine-precision confirmation)

- The restart-density charge upgrade went live (machine-precision dq under FD; the printed-pops
  4-decimal floor bypassed). With it, direct/(µ·q) tracks **√CN(mean-rc)** at ratio 1.07–1.21
  over R = 2..4 — every alternative (plain CN in either rc convention, single-erf refits) varies
  by FACTORS. The µ CN-dependence is √CN-structured, the same family as repulsion's Eq. 56;
  the SI's Eq. 84 prints plain CN_A — the third confirmed paper-vs-binary divergence (Eq. 55
  factor 2; Eq. 84 √; the onsite third-order kernel's degree-2 form).
- OPEN, quantified: amplitude ≈ 1.49× the file's k1 slot value, and a smooth 6% residual drift
  (slightly different counting radii suspected). Separation experiment designed: stretch H
  against DIFFERENT partners (HF, HBr) — different rc, same H constants — isolating rc from
  amplitude cleanly.

### 2026-07-16 (twentieth push) — the DIVERGENCE REGISTRY with mathematical verdicts

- **`PAPER-VS-BINARY-DIVERGENCES.md`** (per Richard): every confirmed paper-vs-binary split,
  with evidence, a mathematical analysis of which side is "correct" in what sense, and the
  implementation verdict. D1 (Eq. 55 factor 2): both means are valid parameterization choices,
  neither derivable at kexp=1.5 — binary correct-by-construction, paper likely dropped the 2.
  D2 (onsite 3rd-order): NEITHER side satisfies strict DFTB3 derivative-consistency (which
  would give a U-INDEPENDENT equal-U onsite); the binary's degree-2 form equals the offsite τ
  structure with the distance factor regularized — our port may offer a labelled
  theory-consistent option. D3 (Eq. 84): family consistency (√CN in Eqs. 56 and 28) favors the
  binary's √CN — a dropped radical sign; amplitude ~1.49× still open.
- Multi-partner round, first attempt: HF's CN window closes below R≈1.9 (F's small radius) —
  probe distances corrected; Br BLOCKED on its reference occupations (SI table ends at Ar) —
  the FD refocc measurement for Z>18 goes first. Both queued.

### 2026-07-16 (twenty-first push) — Br refocc measured (sum exactly 7.00000); multi-partner µ-CN data

- **Br's fractional reference occupations FD-measured** — the first element beyond the SI's
  table: s 1.89888, p 5.02125, d 0.07988, sum 7.00000 EXACT (the electron-count identity is the
  built-in check). Banked; the technique now covers any element the file parameterizes.
- **Multi-partner stretches (corrected windows)**: per-pair the √CN shape holds (HBr's k_eff
  constant to 0.3% over its scan) but the implied pair radii (HF 0.99, HCl 1.47, HBr 1.82) match
  NO decoded radius set, and k_eff varies by partner (1.22/1.29/1.50) — the µ-CN either carries
  its own hardcoded radius table (like refocc and avg_cn) or a second CN-like term overlaps the
  signal. Dataset banked in derived-constants; next discriminator designed: dense per-pair
  R-scans → exact (k_eff, rc) pairs → joint per-element radius solve across many partners.

### 2026-07-16 (twenty-second push) — the µ-CN's OWN TABLE established; I refocc measured

- **Iodine's references FD-measured** (s 1.92602 p 4.98008 d 0.09390 — sum 7.00000 exact),
  joining Br beyond the SI table.
- **Dense scans (4 pairs × 6 points, restart-precision charges)**: √CN explains 99%+ per pair
  (HBr worst-resid 9e-4). The fitted pair radii (HF 0.983, HCl 1.471, HBr 1.824, HI 2.119)
  exclude EVERY decoded radius set, and fixed-radii/free-steepness refits need per-pair
  steepness (−2.39..−1.44) — no shared convention survives. **Conclusion: the µ-CN carries its
  own hardcoded per-element geometry — the FOURTH hardcoded table** (after refocc, avg_cn, the
  onsite functions). k_eff per pair 1.23–1.49 (~5–10% relaxation contamination; underlying
  ~1.3–1.4× the k1 slot).
- **Closure designed — triangulation**: perturbing the PARTNER's k1CN in mixed X–Y pairs
  measures rc(X,Y) directly; enough pairs over-determine the whole element table. Raw scans in
  `data/mu-cn-scans.json`.

### 2026-07-16 (twenty-third push) — the µ-CN is a PAIR function; ES1 closes 4× tighter

- **Both atoms of a pair share the same coupling function**: the Cl-side measurement on HCl fits
  (keff 1.2484, rc 1.4742) vs the H-side's (1.2929, 1.4708) — the µ-CN count is a PAIR property
  (one bond, counted from both ends), collapsing the "per-element table" into pair radii that
  per-element values must combine to.
- **ES1 forward model with the measured pair C(R)**: worst |d| = 2.4e-3 at R=2 (was 1.0e-2 —
  4× tighter), mid/long range 1–4e-4. Remaining short-range suspects: relaxation haze in keff,
  the f-linearization at |q|~0.33, offsite screening. The first-order term is a whisker from
  closed on HCl.

### 2026-07-16 (twenty-fourth push) — Richard's sequence: ES1/ES23 rounds, GRADIENTS exposed, Hamiltonian anchored

- **ES1 sliver**: the small-q f-refinement hit a physical wall — the long-R extraction regime is
  corrupted by the soft-charge mode (at R=10 HCl ionizes fully, q_H = 1.00000, and f(1) = 1.0161
  REPRODUCES the atomic-ion sample: technique validated, regime unusable). ES1 stands at 2.4e-3
  worst / 1–4e-4 mid+long; the refinement path is documented.
- **ES2+3**: the kU CN-couplings measured on HCl (√CN-compatible, rc 1.46/1.60) — but at ~5% of
  the µ-CN magnitudes they cannot carry the short-range residual; the offsite third-order/
  screening structure is the remaining owner (next experiments named).
- **GRADIENTS EXPOSED (two channels)**: (1) the oracle's `-grad` numerical forces validated
  against FD of its own total (~1e-6 agreement — the reference channel for G4); (2) the PORT'S
  FIRST ANALYTIC GRADIENT — repulsion forces with the full CN chain rule — gated at 3.9e-11
  against internal FD (`repulsion.py --grad-test`). One defaults-mismatch bug caught by the gate
  itself (mean_rc convention in the FD comparison).
- **HAMILTONIAN ANCHORED**: H-atom FD gives dε/dL2[0] = −1.00000 Eh/unit EXACT (and per-electron
  in the energy) — **L2 = the shell levels, minus-sign convention**; L3/L4/L1[7] exactly inert
  for lone atoms (molecular couplings, as the map said). The eigenvalue-gate program has its
  first measured anchor: ε(H) = −L2[0] + exchange/spin potentials (−0.2162 − 0.1918 Eh).

### 2026-07-16 (twenty-fifth push) — the atom Fock diagonal DECODED (O-atom probe round)

One FD round on the oxygen atom returned exact structure for every diagonal contribution:
- **L2 = levels at exactly −1.000 per shell**, both spins, strictly diagonal (multi-shell ✓).
- **The first-order µ enters the Fock** (deps/dµ_l = −1.006 for every orbital of shell l) —
  SI Eqs. 87/88 observed live.
- **The exchange potential is occupation-scaled**: occupied orbitals respond to L5[l] at −0.122,
  virtuals at 0.000 EXACTLY — V_ex ∝ L5[l]·n_orbital with a shell-independent 0.122 factor
  (origin open: one more probe class).
- **The spin potential has cross-shell structure**: sign flips with channel; the s shell
  responds (±0.076) despite zero s-magnetization — V_spin(l) = Σ_l' W_ll'·m_l' with hardcoded
  shell-pair weights scaled by L1[9] (w·m_p: s 0.076, p 0.179 at m_p = 2).
- The atom-eigenvalue MODEL is now assembled up to two constants (the 0.122 exchange factor,
  the W shell-pair weights) — the eigenvalue GATE (predict all 8 O eigenvalues + H's 2 from
  file values) is 2–3 probes away.

### 2026-07-16 (twenty-sixth push) — THE MASTER ELEMENT FUNCTION: exchange and 2nd-order damping are ONE law

- **E_x(atom) = −c_x·Σ_l L5[l]·Σ n²** (self-exchange) closes the printed exchange energy
  (O to 4 digits; homogeneity in L5 verified) AND explains the occupied-only eigenvalue
  response (−2c_x) simultaneously.
- **c_x follows the s-rule's exact pattern** (per-period linear, slopes halving period 2→3) —
  and **c_x = s/9.59 for ALL NINE measured elements including H across three periods**: one
  hardcoded per-period-linear MASTER FUNCTION g(el) underlies both the second-order damping
  (s = A·g) and the exchange self-energy (c_x = B·g), B/A = 0.10427. Two of the binary's
  hardcoded element functions were one function all along.
- **E_spin = −½·L1[9]·Σ w̃_ll'·m_l·m_l'** closes oxygen EXACTLY (−0.0702 vs printed −0.070021);
  H gives w̃_ss = 0.2005, N implies w̃_pp = 0.0797 vs O's 0.0895 — the shell-pair weights carry
  remaining element structure (next charting target).
- The atom-eigenvalue gate is now one W-matrix away.

### 2026-07-16 (twenty-seventh push) — the atom-eigenvalue gate: first run, spin EXACT, one term to fix

- The first assembled gate (H + O, all 10 eigenvalues from file values + decoded structures):
  **the spin physics is exactly right** — H's α/β splitting predicted 7.64 eV vs oracle 7.64
  (both channels miss by the IDENTICAL offset), O's channel splittings likewise consistent.
- The misses are clean per-shell constants (H +3.24 eV both channels; O s +13.34, p +4.5–4.8) —
  the signature of a sign/bookkeeping slip in ONE term (prime suspect: the µ-potential F1's
  sign chain through q = refocc − pop; second: the Vx occupied-branch magnitude vs the −2c_x
  measurement). The gate machinery itself works; the fix is localized.
- NOTE for the fix round: verify F1 = −µ_l·f − (Σµq)·f′ sign-by-sign against dε/dµ = −1.006,
  and re-derive Vx from E_x = −c_xΣL5·n² (∂/∂n_iσ = −2c_x·L5·n_iσ ✓ matches the FD) — then the
  shell offsets should collapse.

### 2026-07-16 (twenty-eighth push) — the ACP was the missing term; layout decoded; gate closure at 99.7%

- The eigenvalue-gate misses identified themselves: the ORIGINAL O-atom sweep had already shown
  L9 (the ACPs) moving the electronic energy — atomic correction potentials are exactly
  per-shell constants for an atom, the misses' signature.
- **L9 layout decoded**: [c_s, c_p, c_d, c_f | ζ_s, ζ_p, ζ_d, ζ_f] — coefficients respond
  linearly and shell-diagonally (O: ∂ε_s/∂c_s = +0.9236, ∂ε_p/∂c_p = +0.8773 Eh/unit; d/f
  slots inert without d/f shells), exponents nonlinearly.
- **The closure check**: κ_l·c_l accounts for the gate misses at 99.7% (s) and 98.4% (p) —
  the ACP is confirmed as the one missing term. κ_l is the projector matrix-element structure:
  measurable per element by FD, or computable analytically with our gated integral engine
  (Gaussian projector overlaps — same machinery as Term 4). Next session: κ in place → the
  atom-eigenvalue gate PASSES → extend across elements → diatomics.

### 2026-07-16 (twenty-ninth push) — THE ATOM-EIGENVALUE GATE PASSES; the ACP is pure integrals

- **κ computed analytically = the measured responses to FOUR DIGITS** (O s 0.9236 = 0.9236,
  p 0.8773 = 0.8773; H 0.9105 → V_acp = −0.1190 = the measured miss EXACTLY): the ACP diagonal
  is V_acp(l) = c_l·|⟨AO_l|g(ζ_l)⟩|² — file values + our gated integral engine, NO fitting.
- **THE GATE**: hydrogen — both eigenvalues at 0.000 eV (EXACT: level + µ-potential + exchange
  + spin + ACP all correct simultaneously); oxygen — six of eight orbitals within 0.03–0.07 eV,
  worst 0.30 eV (the β-occupied p; residuals at the size of the FD-measured W-weights and the
  f-slope approximations — refinement targets named).
- The Fock diagonal of the method is DECODED AND GATED at the atom level. Road: refine the O
  residuals (exact W-weights, f-slope) → sweep the gate across elements → diatomics
  (off-diagonal EHT × diatomic frame × our machine-exact overlap).

### 2026-07-16 (thirtieth push) — the gate SWEEPS: five elements, H exact, error budget fully attributed

- Spin weights measured for C/N/F (one probe each; smooth element trend — own law to chart).
- **THE FIVE-ELEMENT GATE**: H 0.000 eV (exact), N 0.160, F 0.216, O 0.299, C 0.623 worst —
  from file values + decoded structures + the ANALYTIC ACP, zero fitted quantities.
- **The error budget attributes itself**: residuals track the one knowingly-unmodeled term —
  the C0 (fourth-order onsite) potential — largest exactly where C0 is largest (C 0.0091 →
  0.62 eV; N 0.0016 → 0.16; O 0.0001 → its 0.30 sits instead in the β-channel exchange fine
  structure). The two named refinements (C0(q) closed form; β-channel split) are the queued
  quartic charting and one fingerprint probe respectively.

### 2026-07-16 (thirty-first push) — the OFF-DIAGONAL channel opens: H2 Fock inversion

- **The instrument**: H2 has one orbital per atom, so the printed eigenvalue pair + our
  machine-exact S12 INVERT to the complete 2×2 Fock at every distance — F11(R) and F12(R)
  measured directly across R = 1..5 (data/h2-fock.json). The off-diagonal EHT rule, the last
  major undecoded structure, is now a measurable curve rather than an equation to trust.
- First readings: F12/S12 drifts −0.53..−0.79 (not a constant Hückel prefactor → the
  Hamiltonian-basis scaling and/or a distance polynomial participate); F11's long-range limit
  sits 0.141 Eh above the atom anchor, with the shift structure accounted by the
  density-dependent terms (delocalization halves the self-exchange; the closed shell turns the
  spin term off).
- **Decomposition designed**: perturb L5 and re-invert → the exchange part of F12 (L5-linear)
  subtracts out, leaving pure H0_12(R) to chart against K·h·S̃ candidates with scaled bases.

### 2026-07-16 (thirty-second push) — the Hückel prefactor MEASURED: H0 = K·S within 0.6%

- The exchange-subtraction decomposition ran: X12 = L5·(dF12/dL5) is 32–47% of the raw
  off-diagonal, and the remainder obeys **H0_12 = −0.400·S12 to ±0.6% across R = 2..4** (the
  raw ratio drifted 50%). The Hückel rule is K·S to first order, prefactor measured.
- Edges (R=1: −0.360; R=5: −0.422) carry the remaining distance-polynomial / basis-scaling
  structure; the ζ-scaled candidate overlap was slightly WORSE than the primary S — the primary
  basis carries the off-diagonal at this order.
- Next: express −0.400 as K·(h1+h2)/2 with the decoded diagonals → locate K among the remaining
  globals; heteronuclear pairs (HF/HCl Fock inversions on the s-block) discriminate the
  averaging rule; then the edge polynomial.

### 2026-07-16 (thirty-third push) — the off-diagonal COMPOSITION table; G1[2] numerology executed

- The K = G1[2] hypothesis (a 0.5% numerical match!) was killed properly: dF12/dG1[2] = 0.00000
  exactly. Measurement over numerology, once more.
- **K = 2.26 measured directly** (the element-L2 perturbation shifts BOTH centers: dF12 = −K·S).
- **The composition table** (Fock-inversion FDs, R = 2.5): the level enters h̄ at weight ~1;
  µ at ~0.45; L3 inert; **L4 is an off-diagonal-ONLY parameter** (dF11 = 0); **L1[7] is purely
  off-diagonal and the strongest mover** (−1.17; the SI's "decouple overlap from Hamiltonian"
  scaling, now localized); **the ACP has a two-center off-diagonal part** (+0.565 =
  ⟨AO1|g⟩⟨g|AO2⟩ — analytic with our engine, like the diagonal κ).
- Structure: F12 = EHT-core (K, levels/µ, L1[7]/L4) + ACP_12 (analytic) + X12 (exchange,
  separable). Next session: compute ACP_12 analytically → subtract → decode the EHT core's
  exact form → the diatomic gate.

### 2026-07-16 (thirty-fourth push) — the EHT core extracted; the Hamiltonian basis is DIFFUSE

- **The F12 decomposition is complete and self-validating**: exchange (slot-linear FD) +
  two-center ACP (analytic — matches its FD response EXACTLY, −0.0738 vs −0.0739) + the naked
  EHT core, extracted across R = 1..5.
- **The tail identifies the Hamiltonian basis**: EHT12/S̃12 with exponents × 0.6 is constant to
  ±1% for R ≥ 2.5 (plateau −0.162) — the Ham basis is a DIFFUSE rescaling of the primary basis
  (consistent with the SI's "scaled exponents… decouple overlap from Hamiltonian"). Short range
  needs one more ingredient. The 0.6 ≈ L4²/2 numerical closeness is FLAGGED AS HYPOTHESIS ONLY
  pending an FD test (the G1[2] lesson).
- Next session: fine-scan the tail scale; FD-test the L4/L1[7] links to the scale; identify the
  short-range ingredient; heteronuclear inversions; then the diatomic gate.

### 2026-07-16 (thirty-fifth push) — the three probes ran: L4 IS the Hamiltonian-basis scale (FD-CONFIRMED)

- **The push-34 one-off is now a committed instrument** (`prototype/eht_h2.py`, labeled data in
  `h2-eht-probes.json`): selftest reproduces the banked R=2.5 Fock row to all digits, the
  analytic both-center ACP matches FD on BOTH matrix elements (5e-4 on the tail), L5 curvature
  1.3e-5. The old unlabeled `h2-h0.json` col4 retro-identified (raw L5 slope, 15-digit match).
- **PROBE 1** (fine-scan): k* = 0.602, plateau −0.1641, tail spread 1.79% (R = 2.5–6).
- **PROBE 2** (the FD test the G1[2] lesson demands): perturbing L4(H) by +0.10 moved k* to
  0.702 vs the hypothesis prediction (L4+δ)²/2 = 0.713 — **inside the pre-declared band. The
  0.6 ≈ L4²/2 numerical closeness is now a CONFIRMED law**: the Hamiltonian basis is the primary
  basis with exponents × L4²/2. Whether the rule is per-AO or per-pair awaits heteronuclear.
- **L1[7] arm**: k* pinned (0.608), plateau ×1.040 for +10% — an AMPLITUDE channel, not a scale.
  SET-to-zero: its content is 26–36% of the core, shaped like the diffuse S̃ itself (subtracting
  it DEGRADES the plateau 1.79% → 6.72%): part of K·h̄, keep it inside. Slope −1.1703 at base
  (linear δ 0.02–0.05, F11 dead), superlinear far out. A prose misread this session (0.8142 —
  that is L1[6]=kU) voided the arm's original k*-band before the run; the discriminator actually
  used (k* moved vs pinned) is value-independent and the correction is recorded in the docstring.
- **ACP relaxation channel isolated**: the c_s FD response exceeds the bare analytic element by
  +0.6% (R=2.0) → +4.4% (R=0.9), zero on the tail — SCF density relaxation; the decomposition
  subtracts the BARE element by construction.
- **PROBE 3** (short-range ingredient): pre-declared bar NOT met — honest abstention. CN is the
  best 1-parameter suppressor (rms 0.0032, ~50% suppression at CN≈0.8) but only 1.17× under
  √CN, and one homonuclear stretch cannot de-collinearize the candidates. Next instrument:
  heteronuclear Fock inversions (HF/HCl s-block) — they separate the candidates AND pin the
  level-averaging rule and the L4 pair-vs-AO question in the same runs; then the diatomic gate.
- Housekeeping: WSL home had lost ~/.basisq and ~/.eeq (oracle failed on startup) — restored
  from /opt/gxtb-v1; the pristine parameter file verified byte-identical to the release asset.

### 2026-07-17 (thirty-sixth push) — the FULL Fock reconstructed; F2 breaks the naive law open; the metric instrument

- **fock_recon.py, gated**: F = S·C·diag(eps)·C^T·S from the restart MO matrix + printed
  eigenvalues + our overlap. Ortho 4.4e-16 (H2)/2.4e-9 (F2), density exact, ties to the 2×2
  inversion at 2.3e-9, and the EXACT π zeros pin the oracle's p ordering = (x,y,z) — something
  the Mulliken gate was provably blind to. Scope measured: the eps print is a window of
  n_occ+2, so H2/HF/F2/H2O invert completely; HCl does not (d shell). Housekeeping: the
  density gate cannot see equal-occupancy column mispairing, but F can — checked at F2 R=6
  (diag ± cross reproduces the printed σ pair) before trusting anything.
- **F2 first, HF second, by design**: homonuclear q=0 kills the polar ES off-diagonal exactly
  as in H2, and fluorine's two shells (L4_s 0.7246 vs L4_p 1.1269) carry the per-AO-vs-per-pair
  question inside one molecule. HF's residual will later BE the ES measurement.
- **f2_stretch.py verdicts (pre-declared)**: per-AO transfer of the L4²/2 law REFUTED (nothing
  flattens, spreads 81–147%; 2-D argmin far from prediction); onsite s-pσ null PASS (3.8%);
  center-swap 2.1e-12; ACP matrix (V diag(c) V^T, all projectors, zero case-work) tail-validated
  2.9e-3.
- **The pσ anomaly and its owner**: EHT(pz,pz') keeps a Coulomb-tailed +0.065 Eh at R=6 — a real
  3.9 eV σ split sitting in the printed eigenvalues, with a FALLING decay rate no Gaussian
  overlap can produce. Seventeen slots FD-excluded (U, µ, offsite-ES, multipole, exchange
  globals, kpen, levels, L4_s, L3); **L4_p owns it** (slope −0.719), L1[7] participates.
- **hmetric.py — ask the binary for its own metric**: S_eff = −(dF/dL2)/K. The LEVEL channel
  rides the PRIMARY overlap (H2 Seff/S = 0.99–1.05) — so the push-34/35 "diffuse Hamiltonian
  basis" is the effective shape of a SUM, not a basis replacement; the L4²/2 FD fact stands,
  reinterpreted as the L4-owned SECOND channel's signature. **F2's s-pσ element CLOSES on the
  metric** (flat ±8% full range, ±3% mid) — the first multi-shell element fully explained.
  Structure now explicit: EHT = level·S_eff + channel(L4, L1[7]); on H2 the channel is the old
  26–36% L1[7] content, on F2-pσ it is the Coulomb-tailed object.
- Next: chart channel_el(R) = EHT − T_lvl·S_eff across H2/F2 (+HF) and decode its closed form —
  ONE undecoded off-diagonal object remains. Then the level weights (w_s ≠ w_p by 6–13%), then
  HF's polar ES, then the diatomic gate.

### 2026-07-17 (thirty-seventh push) — the EULER HOMOGENEITY LAW; the last object cornered and mapped

- **offdiag_atlas.py**: the perturbation map extended to off-diagonal Fock elements —
  rho-normalized slot responses (K and metric cancel). Levels partition EXACTLY by shell pair;
  on (s,pσ) the two level weights SUM TO 1.0000 at every R (0.478+0.522…) — a true weighted
  mean, p slightly heavy. L3/L1[8] are h̄ members with CN(R)-shaped weights — the short-range
  h̄ structure probe 3 was hunting. L4/L1[7] grow with R on every element: channel knobs.
- **THE LAW (closure v2, f2_tail.py)**: the core is DEGREE-1 HOMOGENEOUS in (level, µ) —
  Euler over (L2, L7) closes H2 ss to ±0.006 Eh and F2 ss/spσ/pπ to ≤0.004/0.031/0.013 across
  13 distances. THREE of four element types are now fully explained with no model and no fit.
- **The one survivor** (F2 pσ-pσ): +0.093 Eh at R=2 → +0.064 at 6 → REAL zero crossing
  (σ split smooth through it) → −0.041 at 10 with |obj|·R → 0.41 Eh·Bohr. Its complete
  parameter map: L4_p owns the shape; THREE NEW GLOBALS G1[6]/G1[7]/G2[8] (σ-only:
  spz/ss/π dead to 3e-4 at R=6; G1[7] R-flat) plus L1[7]; ~35 other slots FD-dead. The zero
  crossing + Coulomb asymptote fit a TWO-PART form (overlap-shaped + minus 0.41·γ) — flagged
  hypothesis; 0.41 ≈ 6·c_x(F) is numerology, unbelieved per the G1[2] rule.
- **SI cross-read (Sec 1.3/1.7/1.8)**: Eq 64 H0 = K̄·H̄·Π·S̃sc — diatomic-frame σ/π/δ-scaled
  overlap (Eq 31: harmonic element pairs), shell polynomial Π growing linearly in R (Eq 67),
  shell-wise ζ exponent scales (= L4, our ²/2 convention) and element k̃ charge-adaptations
  (= the L1[7] family) — ALL H̄-multiplied, so Euler kills every Eq-64 form: the object lives
  OUTSIDE Eq 64. Sec 1.8: the g-xTB ACP is NON-LOCAL ONLY — our projector decode is the
  COMPLETE ACP (the ECP-local hypothesis executed by reading, not by fitting).
- F2 tail extended to R=10 (gap-gated, all healthy ≥2.46 eV); hbar FDs completed at all 13 R.
- Next: L8[0..3] as k^diat σ/π/δ candidates (probe at SHORT R — the nullsweep's R=6 blinds
  S̃-shaped responses); locate Eq 67's k^shp element slot; fit the two-part object form; the
  H2-absence constraint on its amplitude (p-specific? drho-carried?); then HF's polar ES and
  the diatomic gate.
