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
