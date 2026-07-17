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

### 2026-07-17 (thirty-eighth push) — the object has a CLOSED FORM; the EHT globals are NAMED

- **atlas_globals.py** (all 20 globals + the L8 quartet at short R — the R=6 null sweep was
  blind to overlap-shaped responses): **G1[0] = k^W_s and G1[1] = k^W_p** (Wolfsberg
  prefactors: rho-constant, exact selectivity — G1[1] dead on ss and on ALL of H2; equal on
  σ and π as Eq. 64 demands). **L8[0]/L8[1]/L8[2] = k^diat σ/π/δ** (Eq. 31): the π slot moves
  pxpx ONLY with exact zeros elsewhere, and hydrogen's π slot is stored as 0.0 — H cannot
  π-bond. **L1[7] = k^shp** (the shell-polynomial element amplitude): its rho is LINEAR in R
  on both molecules — the Π = 1 + k·R derivative signature; 'hbasis' retired. The k^shp,l
  angular-momentum globals are NOT in the file (no global shows rho ∝ R on ss) — a fifth
  hidden-table candidate. Old energy-side labels for G1[0]/L1[7] were these terms'
  through-density footprints; superseded.
- **The pzpz object has a measured closed form** (prototype/object_fit.py, reproducible):
  obj(R) = 0.42·erf(0.23·R)/R − 0.048·S̃(0.635), rms 4.2e-4 Eh on the positive branch and
  4.9e-4 on |obj| over all 13 points — both under the pre-declared 2e-3 bar. The R≥7 sign
  flip is an ATTRIBUTION artifact (the crossing-capable fit chose no crossing; my earlier
  'real crossing' call is corrected). **a = 0.23 = G2[8] (0.2347)**, corroborated by FD shape
  at two distances — G2[8] is the penetration range global. **c = 0.42 ≈ 3·|Δρ⁰_p(F)| = 0.432**
  (3%): FLAGGED numerology that survives its first falsification — REFOCC(H) = 1.0 exactly,
  so Δρ⁰(H) = 0 predicts the object's measured H2-absence. Interpretation (flagged): the
  reference-density PENETRATION correction the SI admits is 'not fully captured' — σ-selective
  because pσ points along the axis, outside Eq. 64, with G1[6]/G1[7] as amplitude knobs
  (entry mode open — G1[7]'s R-flat response is not explained by the c-pathway).
- Next: chart c across elements (the 3·Δρ⁰ test needs a second p element — CO/N2 partial
  inversions or HF's s-p σ), close G1[6]/G1[7]'s entry mode by SET-scans, the H0 forward
  model (all pieces now named: K^W·H̄·Π·S̃sc + object), then HF's polar ES and the diatomic
  gate.

### 2026-07-17 (thirty-ninth push) — the µ-channel LAW; the first forward gate FAILS honestly

- **Eq. 64 read EXACTLY from the PDF** (layout mode): H0 = [(k^W_A+k^W_B)/2]·[(H_A+H_B)/2]·
  Π·S̃sc off-diagonal; same-atom blocks purely diagonal (orthonormality kills onsite
  off-diagonals — consistent with our measured 3.8% onsite null being NON-H0).
- **A decoded law (h0_forward.py part 1)**: the ES1 off-diagonal channel = −µ̄·S_primary with
  coefficient −1.0000 MEASURED to 3–4 digits at long R on every element type — the Mulliken
  ½S(v_A+v_B) form at exactly unit weight. The old drifting µ-rho is fully explained: µ rides
  the PRIMARY overlap while the levels ride Π·S̃sc; the short-range deviation is the µ-CN
  (4th-table) + f(q) structure.
- **The polynomial slopes fall out of banked data**: 1/(rho_L17/R) is linear in R → b(H2) =
  +0.033; fluorine's k^shp = −0.0040 makes Π(F2) ≈ 1 — F2's Seff fits must run b≈0, and the
  free-fit b = +0.14 was compensating an S̃-shape error (caught).
- **THE FIRST FORWARD GATE: FAILED** — pre-declared |diff| ≤ 0.01 Eh at 13 points; worst
  0.028. Causes named: (a) the (C, b, k) Seff fits are degenerate, so amplitude checks
  against k^W·k^diat are not robust from shape fits; (b) under the honest-b constraint H2's
  metric misses ~1% systematically — the Ham basis's k̃-adapted CONTRACTION coefficients
  (SI Eq. 28 tilde set) are unmodeled and live exactly where CN varies; (c) C_µ extrapolated
  below R=1.4. One tantalizing number, NOT claimed: F2-ss constrained-fit C = 0.6680 vs the
  k^W_s·k^diat_σ(F) prediction 0.6708 (0.4%) — the degeneracy forbids the claim.
- Next: measure S̃sc SHAPE-FREE as Seff/Π(known), match basis constructions (density-adapted
  vs k̃-adapted vs scaled exponents) against the measured curve, locate the k̃ tilde slots,
  then reassemble and re-gate. The µ-law and the object's closed form slot straight into the
  eventual assembly unchanged.

### 2026-07-17 (fortieth push) — THE FIRST FORWARD GATE PASSES: H2's off-diagonal predicted from the file

- **The degeneracy breaker**: Eq. 28 (read clean from the PDF — no screenshots needed; layout
  mode extracts the SI equations legibly) says the Ham basis adapts k̃0/k̃2/k̃3, and at q = 0
  only k̃2·√CN survives. Adding that kb dimension to the metric match (M = Seff/Π vs
  S̃(k, kb)) collapses v1's degenerate fits into a SHARP unique minimum: rms 9.9e-5 = 0.01%
  of range (bar was 0.5%), rms growing 50–100× within ±0.05 in k or ±0.10 in kb.
- **The L4 law corrected**: k = 1.110 ≈ L4(H) = 1.094 DIRECT. The old k = L4²/2 was an
  artifact of fitting the level+µ mixture; the probe-2 single-δ FD proved L4 drives the scale
  but could not separate the functional forms. Correction recorded in the anchor.
- **kb = +0.17** (vs the density basis's 0.2272) — the Ham basis's own √CN coefficient,
  measured; its slot not yet located (G1[4]/L8[3]/L8[7] candidates, FD-vs-analytic pending).
- **The amplitude assembles**: C·2.26 = 2.082 vs k^W_s·k^diat_σ(H) = 2.091 — 0.4% (the same
  0.4% F2-ss showed independently). Eq. 64's prefactor is now quantitative from named slots.
- **THE GATE**: with C_µ(R) measured at 7 points (new short-R FDs; exactly −1.000 at R ≥ 4)
  and Π's b = +0.0315 from the polynomial's own response, the forward assembly
  H0 = k^W·k^diat·(−L2)·Π·S̃(k, kb) + µ·C_µ·S predicts H2's measured EHT at ALL 13 points,
  worst |diff| 0.0094 Eh — **PASSED at the pre-declared 0.01 bar** (stretch 0.006 not met;
  residual structure named: the dropped L3·CN term, the 0.4% amplitude gap, sparse C_µ
  between 2.5 and 4).
- Next: locate the k̃2 slot (FD vs analytic ∂S̃/∂kb shapes); the same decomposition on F2's
  shell pairs (Π ≈ 1 there); the σ elements with the penetration object added; HF (the k̃0/k̃3
  q-channels wake up on a polar molecule); then the full diatomic gate and SCF assembly.

### 2026-07-18 (forty-first push, overnight autonomous run) — L8[3] = the k̃2 knob; F2-ss forward PASSES; the polar wall found

- **kb-slot hunt (kb_slot.py)**: L8[3] IS the Ham-basis √CN adaptation knob — its response
  tracks the analytic ∂F/∂kb at constant ratio (±5%) where every alternative swings by
  factors; L8[2]'s numerology executed (exact zeros). kb = 0.4223·L8[3] (linear-through-zero
  consistent: predicts 0.179 vs fitted 0.17–0.18); the 0.4223's globality flagged.
- **F2 forward v3** (channels now MEASURED at all 10 points, f2_channels.py): **ss PASSES**
  (0.0092; k_s = 0.730 vs L4_s = 0.7246 — L4-direct on a second shell at 0.7%; amplitude
  1.0% from k^W·k^diat). The p-elements keep structural short-range gaps (spz worst 0.053):
  k_p fits 7% above L4_p with amplitude 12–17% high — the π/p channel carries something
  unidentified; recorded, not tuned.
- **HF (hf_stretch.py + hf_analysis.py): the polar wall.** On a polar molecule every FD
  response carries first-order density relaxation (homonuclear symmetry had pinned q = 0):
  the Eq-64 weight tests fail with R-drift (µ H/F ratio 1.7–2.0 vs the 1.004 f(q)
  prediction). The FD-decomposition's clean domain is homonuclear; the polar layer must be
  closed SELF-CONSISTENTLY. The forward residual (+0.07 Eh at short R, not ∝ S·q) is banked
  as the SCF stage's target. Pivoting to the SCF assembly (the roadmap's endpoint anyway):
  read Sec 1.9.2/1.10.2/1.15.2 Fock forms, assemble H2 first, gate on converged observables.

### 2026-07-18 (forty-second push, overnight) — THE OBJECT WAS THE EXCHANGE; D4; the kernel measured

- **D4 registered** (mfx_chart.py): the binary's onsite Ex is EXACTLY linear in L5 —
  dEx/dL5 = −0.049268 constant to SIX DIGITS over L5 = 0.5–8, and it equals the master
  function c_x = s/9.59 — refuting Eq. 149's denominator-U printing. Offsite Ex SATURATES
  in L5 (numerator AND denominator). The onsite channel is the separate Sec-1.16 correction,
  hardcoded scale: the all-20 global sweep gives EXACT ZEROS on the atom's Ex (α hardcoded).
- **The exchange kernel MEASURED** (ex_extract.py): H2's printed Ex + restart P + our S
  invert Eq. 151 exactly (the idempotency invariants A = B = −0.2500 at every R confirm the
  Mulliken structure); γ_off(R) = 0.437→0.466 (max at R≈3) →0.377 at R=20 — glacial decay,
  γ·R still rising at 20. The MNKO form-fit reaches 9.9e-4 but is grid-edge degenerate over
  this window: the CURVE is the deliverable. Offsite globals found: G1[4] (strongest),
  G1[6], G1[7], G2[8] — four movers for the SI's four kernel globals.
- **THE OBJECT REINTERPRETED — mystery closed**: the pzpz object IS the offsite-exchange
  Euler remainder. The L5-Euler subtraction is exact only for degree-1 content; the offsite
  kernel saturates, so the subtraction removes the tangent line (+0.027 of +0.077 at R=6)
  and the saturation gap survives — with every signature now explained: σ-only (bond-order
  P² weighting; π and s pairs closed), the kernel's slow tail, the exchange globals as its
  knobs, L4_p through the density's s-pσ hybridization. The "penetration"/3·Δρ⁰ readings are
  SUPERSEDED (recorded); the 0.42·erf(0.23R)/R form survives as the curve's parametrization.
- Consequence for the assembly: the SCF's exchange Fock = Eq. 153 with onsite-linear +
  offsite-kernel-curve pieces; the F2 p-element forward gaps (spz/pzpz/pxpx) are expected to
  be largely THIS remainder — re-gate after the exchange is forward-modeled.

### 2026-07-18 (forty-third push, overnight) — the first SCF: exchange ENERGY exact, the Fock's off-manifold form is the last piece

- **scf_h2.py**: H2 assembled from gated parts (H0 forward + µ-law + analytic ACP + Eq-151
  exchange with the measured kernel). **The exchange ENERGY gates EXACTLY** (d = −0.00000 at
  all three R) and ES1 gates; the occupied eigenvalue lands within 0.008–0.020; **the VIRTUAL
  misses by +0.18–0.33 — gate FAILED through the virtual, honestly.**
- **The diagnosis is clean**: on H2's idempotent manifold the Mulliken exchange energy is
  CONSTANT (the A = B = −0.25 invariants), so no energy check can pin the Fock; my Eq-153
  translation is provably the exact gradient of my energy form (a numerical-gradient Fock
  reproduces it) — the binary's Fock is the gradient of a DIFFERENT off-manifold
  continuation of the same on-manifold energy.
- **The required-Fock curves are extracted** (X_req = reconstruction − H0 − µS − ACP, 12 R):
  X_req_11 → −0.1790 at dissociation = −γ_on/2 EXACTLY (the spin-restricted atom limit);
  my form's limit is wrong (keeps falling). The binary's placement is simpler; fitting its
  form against these curves is the next push's job — then the SCF re-gates, then F2/HF.

### 2026-07-18 (forty-fourth push, overnight close) — the exchange-Fock placement: leading candidate found, floor reached

- The γ∘(P-sandwich) family REFUTED for the Fock placement (unphysical coefficients, wrong
  dissociation limits). The **Mulliken-potential form** −½S∘(v_µ + v_ν), v = γ·(per-spin
  Mulliken populations), matches BOTH H2 elements at a consistent w ≈ 0.47–0.52 for R ≤ 3 —
  its S-suppressed cross element is exactly the required signature — with a small
  slower-falling F12 remainder at long R.
- **The floor**: fits plateau at rms 1.4e-2 because X_required is extracted THROUGH the
  assembled H0+µS+ACP, whose own gaps are ~0.01 Eh. The placement cannot be pinned tighter
  than the pieces beneath it. NEXT (the morning round): refine H0's (k, kb, b) on finer
  grids + C_µ everywhere + the L3 term with the proper internal CN → re-extract X_required
  at the 1e-3 level → pin w and the small term → re-gate the SCF → then F2/HF SCF and the
  energy bookkeeping gate.
- Overnight tally (pushes 41–44): L8[3] = the k̃2 knob; F2-ss forward PASSES; the polar wall
  mapped (FD's clean domain = homonuclear); D4 registered (the exchange kernel's U placement,
  onsite exactly linear at six digits = c_x = s/9.59); the offsite kernel γ(R) measured at
  17 points with its four globals located (α hardcoded); THE PZPZ OBJECT SOLVED (the offsite
  exchange's Euler remainder); the first SCF assembled with the exchange ENERGY gating
  exactly; the Fock placement's leading candidate identified. Every step committed.

### 2026-07-18 (forty-fifth push) — the exchange Fock DECODED to its skeleton: population form, exact atom anchor

- **The manifold-breaking instrument**: H2⁺'s empty β channel makes F^α − F^β = X + spin
  (the spin part measured by its own W-slot FD) — the exchange Fock measured DIRECTLY, no
  forward stacking. The energy functional was already three-manifold-confirmed (H2 exact;
  H2⁺ = −0.125·(γ_on+γ_off) to all digits; H⁻ = −γ_on exactly).
- **The DISCRIMINATOR**: X11(H2⁺)/X11(H2) = 1.01–1.05 — the Mulliken-POPULATION form
  (predicts equal), refuting every P-matrix-linear form (predicts one half).
- **The atom anchor is EXACT**: X^α(H atom) = −0.360142 = −γ_on to five digits — weight 1,
  m^α = 1. And H2⁺'s dissociation limit puts the diagonal's γ_off weight at 0.003 ≈ ZERO:
  the Fock diagonal = −γ_on·(own same-spin population), full stop, plus a small S²-shaped
  short-range correction (ratio ≈ −0.065 (H2⁺) / −0.049 (H2) — the one sub-term left).
- Off-diagonal first reading: X12 ≈ −½S·γ_on·(m1+m2) within 10–18%. Next: the joint lsq
  across (atom, H2⁺, H2) pins the S² term and the off-diagonal exactly → rebuild fock_x →
  re-gate the SCF → F2/HF.

- Addendum (same push): the off-diagonal X12 is NOT S-proportional (a P12-carried piece
  dominates at long R); the 4-term {S, P12}×{γ_on, γ_off} fit reaches 3.1e-4 on the clean
  H2⁺ curve but with degenerate coefficients — measured, not identified. Next instruments:
  joint cross-system fits (H2-neutral + F2 spin probes) or the diagonal-first SCF re-gate.

### 2026-07-18 (forty-sixth push) — the triplet manifold: m-scaling confirmed, a clean triplet law, the algebra task queued

- **h2triplet_fock.py** (the 4th manifold; β empty again, α holds BOTH MOs: m = 1,
  P12 = −S/(1−S²) sign-flipped): the diagonal's m-scaling CONFIRMED (X11 → −γ_on exactly at
  dissociation with m = 1, completing the atom/H2⁺/triplet ladder m = 1, ½, 1), and a clean
  law: the triplet's diagonal extra = −0.154·S·P12, constant ±2% across the stretch.
- The four-manifold joint subset fits (atom + H2⁺ + triplet + singlet, closed-form
  descriptors) reach rms 4–5e-3 with system-structured residuals — the product basis is
  still not the binary's true form. Recorded, not tuned.
- NEXT: the algebra task — derive the Ref-62 four-matrix Fock formula's 2×2 element
  structure symbolically per manifold (the binary's loop conventions leave few candidates)
  and match against the measured curves; then rebuild fock_x, re-gate scf_h2, F2/HF.

### 2026-07-18 (forty-seventh push) — the great elimination: the whole SPS-sandwich space is dead; two new asymptotic laws

- **The six-pairing solve** (all six possible kernel index-pairings in F = −ΣSPS·γ_pair,
  free coefficients, 57 equations over four manifolds): rms 5.2e-2, unphysical coefficients,
  diagonals missed by 0.08–0.13 — the ENTIRE sandwich space is eliminated, every loop
  convention. The Mulliken-density form γ∘(PS+SP)/2 also dies (triplet off-diagonal would
  vanish; measured −0.48).
- **Two new asymptotic laws**: X12(triplet)/S → −γ_on exactly at long R; X12(singlet, R≥4)
  = −γ_on·S·m̄ − 2c_x·γ_off·P12 within 0.6–2% — the 2c_x = 0.0985 factor appears uninvited
  (flagged; c_x is the onsite exchange scale, so a mechanism is plausible). Short-R still
  refuses every closed form tried; recorded, not tuned.
- Next: the sympy symbolic-enumeration instrument over the wider form space (population,
  Mulliken-density, P-carried, and product structures), solved as exact rational identities
  per manifold against the four measured curve families.

- Addendum (same push): **the onsite-kernel Fock part is IDENTIFIED** — F^X(onsite) =
  −½S∘(v+v), v_A = γ_on·m_A^σ — unifying the atom anchor, the m-ladder, and both molecules'
  off-diagonal asymptotics in one form. SCF re-gate with it alone: **R=2.5 passes fully**
  (both eigenvalues ≤0.006); worst virtual miss down 0.33 → 0.10. The remaining misses are
  the measured γ_off remainder, provably two pieces crossing at R≈2.3 (long branch =
  −2c_x·γ_off·P12, flagged). The symbolic session takes the remainder; then the full gate.

### 2026-07-18 (forty-eighth push) — the remainder characterized; the singlet's short-R data ruled contaminated

- Trust correction: the singlet's short-range remainder is FORWARD-CONTAMINATED (it rides
  the H0-forward's ~+0.01 gaps); the clean set is triplet + H2⁺ (spin-subtracted) + the
  singlet tail. On it: H2⁺'s remainder = the −2c_x·γ_off·P12 law ±15% (R ≥ 2.5) plus a
  small short piece; the triplet's is one negative fast-dying piece that resists every
  one/two-factor product tried (all drift >25% — recorded, not tuned).
- The symbolic-enumeration session now has everything it needs: three clean curve families,
  the identified onsite form, the lawful long branch, and the crossing structure. Then the
  full H2 gate → F2/HF → energies → the port.

### 2026-07-18 (forty-ninth push) — two exact laws: the triplet Fock is PURE Mulliken-potential; the (P−m) object

- **rem(triplet) = ρ(R)·S EXACTLY** (rem12 = rem11·s to 3–4 digits everywhere): the triplet's
  entire exchange Fock is the Mulliken-potential form with ONE extra per-atom potential
  w = ρ/2. And w tracks **(P_AA − m_A)** — the onsite density-matrix element minus the
  Mulliken population — at c ≈ −0.080; H2⁺'s short-range OFF-diagonal piece fits the SAME c
  (−0.0802 vs −0.0815). H2⁺'s diagonal short piece still resists (sign conflict) — recorded.
- The hybrid hypothesis (onsite population-gradient + variational γ_off-gradient) REFUTED
  cleanly: parallel-matrix ratio ≈ −1.33 on the triplet (wrong sign), wrong shape on H2⁺.
- The Fock's emerging complete form: −½S∘(v+v), v_A = γ_on·m_A + w_A(P, S) with w built
  from (P_AA − m_A) and P12-carried pieces (coefficients ≈ −0.08 and −2c_x — neither
  identified against known constants yet). One unification step left.

### 2026-07-18 (fiftieth push) — the (PS−I) object; the purity diagnostic; a correction

- **Purity diagnostic**: pure Mulliken-potential ⟹ X12 = s·X11 (symmetric systems). The
  triplet passes EXACTLY; H2⁺ violates 8× at long R → the complete form = −½S∘(v+v) + an
  off-diagonal-only term. **That term's carrier is (PS − I)**: identically zero for the
  triplet (P = S⁻¹), exactly ½ off-diagonal for the paired manifolds at all R — one object
  explains the triplet's purity AND the R-flat H2⁺ branch (≈ −2c_x·γ_off·(PS−I)₁₂, ~7%).
- Triplet potential sharpened: w = −0.166·γ_off·(P−m), constant ±3%; −0.166 ≈ −⅙ (flagged).
- **Correction recorded**: the earlier "same c on H2⁺" claim skipped the Mulliken chain's
  −s factor; corrected, H2⁺'s short pieces are LAWLESS still (and the two off-diagonal
  pieces cancel near R=2.3, so the SCF needs both or neither).
- Next discriminator: HeH⁺ — heteronuclear, one electron: two elements' kernels in one
  clean spin-subtracted system; strong test of every γ-weighting question at once.

- Addendum (same push): **HeH²⁺ is viable and He is calibrated** — Ex(He) = −0.8992 =
  −γ_on(He) with dEx/dL5 = −0.19224 = −2c_x(He) EXACTLY (self-consistent): c_x(He) =
  0.09612, s(He) = 0.9218 — a NEW element datum giving the s-rule its period-1 line
  (a = 0.4492, b = 0.0234). The heteronuclear spin-subtraction campaign opens next session.

### 2026-07-18 (fifty-first push) — heteronuclear inversion bug caught; the UKS restart DECODED

- **Bug caught before it bit**: the symmetric 2×2 eigenvalue inversion assumes F11 = F22 —
  invalid heteronuclear (two eigenvalues cannot determine three unknowns). heh_fock.py's
  first-run X columns are flagged INVALID in the JSON (the giveaway: "F12" = −0.44 between
  decoupled fragments). All homonuclear results stand — symmetry guaranteed their inversion.
- **The UKS restart layout is decoded and verified**: rec0 = P_α‖P_β packed (β all-zero for
  one-electron systems), rec1 = C_α‖C_β column-major — H2⁺'s P11 = 1/(2(1+s)) and the
  bonding/antibonding columns reproduced exactly. The per-spin full reconstruction
  F_σ = S·C_σ·ε_σ·C_σᵀ·S is unlocked; the heteronuclear campaign reruns properly next.
- Also learned: HeH²⁺'s electron sits ~fully on He at every R (m_He = 0.98–1.0) — the
  system is weakly mixing for cross-kernel purposes; alternatives (LiH⁺-type) noted.
- He calibration stands (γ_on(He) = 0.8992, slope-exact); γ_cross(HeH) ≈ +0.06–0.08 mid-R
  (weakly determined, no inversion involved).

### 2026-07-18 (fifty-second push) — the heteronuclear campaign CONFIRMS the skeleton; UKS gate machine-exact

- **The UKS reconstruction gates at machine precision** (H2⁺ ΔF vs the symmetric inversion:
  1.7e-16) — both instruments validated against each other where both apply.
- **HeH²⁺ proper**: the population law holds with ELEMENT-RESOLVED kernels — X_HeHe =
  −γ_on(He)·m_He to ≤0.006 at every R (through m_He > 1, and −0.00004 at dissociation);
  X_HH = −γ_on(H)·m_H exact asymptotically INCLUDING m_H < 0; the off-diagonal
  −½s(v_He+v_H) ≤ 0.02 everywhere. Five manifolds now confirm the exchange-Fock skeleton.
- **New localization**: the residual short-piece lives on the SMALL-population atom (H:
  −0.086 at R=1.2 vs He 30× smaller, same geometry) — leverage no homonuclear system could
  give. The (m−1)-carrier candidate fits here (+0.081, the −0.08-family again) but would
  vanish on the triplet, which has its own extra — the unified short-piece law is still the
  one open item. restart.py + fock_recon.py now carry UKS natively.

### 2026-07-18 (fifty-third push) — the n=2 short-piece law UNIFIED; two refutations

- **Triplet HeH⁺** (2α, heteronuclear, P = S⁻¹ to 1e-8): the short-piece potential is
  ELEMENT-BLIND — wHe = wH to four digits at every R (He's constants are 2.5× H's: every
  per-element weighting refuted) — and the γ_off-weighting is REFUTED 24× (the tiny cross
  kernel predicts −0.012; measured −0.284).
- **THE UNIFICATION**: both triplet systems obey rem11 = −0.45·γ̄_on·(P_AA − m_A) with
  γ̄_on the PAIR-MEAN onsite kernel (HeH⁺ ±0.6%; H2 ±3%); c ≈ −4/9 flagged. The off-diagonal
  keeps the pure Mulliken-potential pattern (rem12 = s·rem11) on both.
- The n=1 manifolds carry a structurally different piece (small-population-localized,
  sign-opposed to (P−m)) — a SEPARATE term, still lawless: the last unknown in the entire
  Hamiltonian is now one n=1 short-range term with four clean datasets bearing on it.
