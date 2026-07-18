# g-xTB Port — Ledger & Decode Status

**Authoritative current state of the parity ledger and the decode work.** Read this instead
of digging through the `TODO.md` journal (which is an accumulating pass-by-pass history). Last
updated 2026-07-18.

The goal: a faithful python/Fortran g-xTB whose gas-phase energy matches the reference binary
(`/opt/gxtb-v1/gxtb`) term-for-term, so our own solvation layer can sit on top of it. The
instrument is `energy_parity.py` — every printed term of the binary's energy decomposition vs
ours, in mEh, over six neutrals (H2, F2, HF, H2O, CH4, NH3).

Run it: `GXTB_V1_EXE=/opt/gxtb-v1/gxtb python3 energy_parity.py` (gpudft conda env, WSL).

---

## The ledger — where every term stands

| term | worst \|d\| (mEh) | status |
|---|---|---|
| **Ex** (exchange) | **0.000** | DECODED + wired, bit-exact on all six (`mfx.py`) |
| **atomic core** | **0.000** | measured constants, exact |
| **nuclear repulsion** | **0.003** | exact (needs the **EEQ** charge, not Mulliken) |
| **dispersion** | **0.299** | wired: stock dftd4 + derived damping (`dispersion.py`); under grade |
| **ES1** (charge SIE) | 0.750 | wired; F2 homonuclear edge is the worst point |
| **ES multipole** (AES) | 0.800 | energy DECODED; the epol self-term kernels are the open piece |
| **ES2+3** | 3.198 | ES2 wired; the **ES3 energy** (water) is the open piece |
| **electronic** (H0) | **321.7** | **the #1 gap.** Now measured on all six; per-element params sourced |
| Espinpol | 0.000 | zero on closed shells |

Per-molecule **electronic** (H0) defect: H2 +47, F2 +16, HF −21, H2O +84, CH4 +295, NH3 +322 mEh.
It grows with p-shells × atom count — dominated by the missing off-diagonal shell polynomial.

---

## Wired / decoded (done)

- **Exchange (Ex)** — the pass-88 range-separated Mulliken exchange, `2·mfx.ex_energy(P,S,mfx.gamma_matrix(...))`. 0.000 mEh everywhere.
- **Restart density** — the ledger takes P straight from the binary's checkpoint (`restart.converged_state`); this is what gave CH4/NH3 full coverage.
- **Dispersion** — `dispersion.py` calls stock dftd4 with the g-xTB damping derived from the parameter file: `a1 = g1[9] = 1.2154627292`, `s8 = g2[9] = 0.304294728`, `a2 = 0` (revD4 eliminates it), `s6 = s9 = 1`. Lands 0.03–0.30 mEh from the binary — under substitution grade. The residual is revD4's Mulliken-charge ζ (SI Eq 166-169), not yet ported. `dftd4>=4.2.0` in `requirements.txt`, present in the gpudft env.
- **electronic measured on all six** — `Tr((H0+A)·P) + Ex`. GE.build already carries a polyatomic H0; with Ex exact everywhere this is a clean H0-defect readout (was hidden as "MISSING").

---

## THE KEY ASSET — the `gp3_` parameter symbol table

The binary keeps **every** method parameter in a named module array (`parcom_mp_gp3_*`), and
the symbol names say which parameter-file column feeds which SI equation. This ends the
column-guessing the port had done by hand and sources the whole per-element layer by name.

Full how-to and the complete map: **`GP3_PARAMETER_SYMBOLS.md`** (dedicated doc) and
`data/derived-constants.json` → `gp3_parameter_symbols_ROSETTA`. Highlights: it sourced all four
per-element H0 pieces (`gp3_lev_`, `gp3_lev_cn_`, `gp3_poly_`, `gp3_rad_`) and confirmed the ES3
Γ table is `gp3_gam3_` (bit-identical to the extracted `rsi`).

---

## Open items, ranked, with exact next steps

### 1. Core Hamiltonian (H0) — the #1 gap (electronic, up to 322 mEh)

**Structure fully decoded** to SI Eq 64-67. GE.build's H0 drops two pieces:
- **CN-dependent shell levels** (Eq 65): `H_lA = L2 − LEVCN·CN`, with `LEVCN = gp3_lev_cn_ = shells[1]`, `CN` = the internal (`es2_energy.coordination`) one. Sign VERIFIED (the flip is 5-9× worse).
- **Shell polynomial** (Eq 66-67): `Π = π_lA·π_lB`, `π = 1 + k^shp·k^(shp,l)·(R/Rcov)`, with `k^shp = gp3_poly_ = l1[7]`, `Rcov = gp3_rad_ = l1[5]`.

**Implemented + measured** (patch: `data/h0_v2_CN_levels_and_poly.patch`): the total error
across the six halved (785→283 mEh), CH4/NH3 improved several-fold — but the polar pair (HF/H2O)
regressed. **Reverted** (GE.build is load-bearing for the SCF; a mixed change with an unsourced
global must not ship).

**What's left:** (a) the 4 `k^(shp,l)` globals — mixed into g1/g2, *not* cleanly separable, must
be extracted not fit; (b) the σ/π `kdiat` treatment (currently a v1 harmonic mean with a crude
switch — the likely cause of the polar regression).

**Recommended route — extract the binary's H0 via the eigensolver:** g-xTB pre-transforms then
calls mkl **DSYEVR** (`0x5a0ee0`, standard symmetric, 3 hits = SCF iterations). Arg `r8` = the
orthogonalized Fock `Xᵀ F X`. Tested this cycle: the *first* call is **not** bare H0+ACP (it
carries the EEQ guess density / a non-Löwdin transform — naive reconstruction missed by ~1 Eh).
So: read A at each of the 3 calls, match eigenvalues to the printed `eps` to find the converged
one, determine the transform (require un-transform to be symmetric), then subtract the electron
terms we already have (ES1/ES2/ES3/Ex/AES). Full note: `derived-constants.json` →
`h0_electronic_gap_ROADMAP.extract_H0_via_DSYEVR`.

### 2. ES3 energy contraction (ES2+3, 3.2 mEh on water)

Everything but the energy is done: V3 is a pure quadratic in the shell charges, the full
interaction tensor is measured on HF + H2O (7-figure), the onsite law is exact
(`τ_onsite = −½·γ2²`), and Γ = `gp3_gam3_` is confirmed. **The blocker:** the naive contraction
`(1/3)·q·V3` and the direct Eq-129b formula both give **−0.18 mEh**, but the true ES3 (printed
ES2+3 minus the re-verified E2 = 0.045023) is **+0.024** — wrong sign. So `set3espot_`'s V3 is
*not* `∂E3/∂q` for the printed energy; the scf assembly (`scf_+41325 = 0x58789d`, where V3 is
consumed) combines it non-trivially. That assembly is the decode target. Probes:
`es3_probe.sh`, `es3_sweep_probe.sh`, `es3_h2o_probe.sh`. Detail: `derived-constants.json` →
`es3_injection_round1/2` + `round3_energy_SHARPENED`.

### 3. AES epol kernels (ES multipole, 0.8 mEh — sub-grade)

Formula decoded: stock `aniso_electro` = `e01(gab3) + e02(gab5) + e11(gab5) + epol(self)`. aes.py
has two bugs — it adds spurious g7/g9 (stock has only gab3/gab5) and omits the epol self-term
`Σ dipKernel·|dip|² + quadKernel·|qp|²`. On the charge-free diatomics the whole defect IS epol.
The kernels are **CN-dependent** (a CN-linear model closes H to <1 µEh). **What's left:** the exact
per-element kernels — no `gp3_` symbol for them (they're in the stock aesData struct, a separate
path); Richard's lead is that the CN they use is the *other* one, not es2. Probe:
`aes_epol_probe.sh`. Detail: `derived-constants.json` → `aes_multipole_energy_DECODED`.

**Do NOT re-apply without the epol term:** removing g7/g9 alone regresses the polar gate (the two
errors currently partially cancel).

### 4. Dispersion ζ (dispersion, 0.3 mEh — sub-grade)

The stock-D4 wiring is under grade; closing the last 0.3 needs revD4's sigmoidal Mulliken-charge
ζ (SI Eq 166-169), which needs the atomic reference polarizabilities α+/α0/α−. Detail:
`derived-constants.json` → `dispersion_revd4_slots_DERIVED`. Probe: `d4_probe.py`.

---

## Methods that worked (Richard's, banked)

- **Read at the converged/last call** — a routine called N times can accumulate or evolve its
  argument across calls (set3espot evolves EEQ→Mulliken over its 9 calls).
- **Reads-only under gdb, never inject** — Fortran block-buffers piped stdout, so any tampering
  that crashes the tail loses the already-computed (buffered) printed line. This defeated every
  array-injection until reads-only.
- **Two coordination numbers** — the binary has an internal covalent CN (ES/H0/MFX/repulsion)
  and a separate basis/dispersion CN; use the right one per term (H0 levels + AES self-term want
  the internal one).
- **Synthetic-input probing** — inject whole/unit/doubled values into a routine's inputs to make
  terms drop out and read the structure (proved V3 quadratic, measured the ES3 tensor).
- **The `gp3_` symbol table** — source any per-element parameter by name (see above).

---

## File index

| what | where |
|---|---|
| the ledger | `energy_parity.py` |
| exchange | `mfx.py` |
| ES2 | `es2_energy.py` |
| AES | `aes.py`; probe `aes_epol_probe.sh` |
| dispersion | `dispersion.py`; probe `d4_probe.py` |
| ES3 probes | `es3_probe.sh`, `es3_sweep_probe.sh`, `es3_h2o_probe.sh` |
| the engine (H0/ES/etc.) | `gxtb_engine.py` |
| H0 v2 patch | `data/h0_v2_CN_levels_and_poly.patch` |
| every decoded constant + all detail | `data/derived-constants.json` |
| the parameter symbol table | `GP3_PARAMETER_SYMBOLS.md` |
