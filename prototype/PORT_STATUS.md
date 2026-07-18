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

**Structure decoded to SI Eq 64-67; three of the four pieces are now SETTLED, and the last is a
diagnosed multi-error tangle.** Detail: `derived-constants.json` → `h0_electronic_gap_ROADMAP`
(entries `v3_*`, `v4_*`, `v5_*`). The v2/v3 implementation lives env-gated in `gxtb_engine.py`
(`GXTB_H0V2=1`, default OFF — v1 ships untouched); patch also at `data/h0_v2_CN_levels_and_poly.patch`.

SETTLED this arc:
- **Shell polynomial (Eq 66-67):** the per-shell factor is the named `gp3_pln_` array (per-element
  `shells[2]`), *not* a hidden g1/g2 global. With it + the CN-levels, **H2 closes +47→+0.3 mEh**, no
  fit. `π = 1 + gp3_poly[Z]·gp3_pln[Z][l]·(R/Rcov)`.
- **CN-dependent levels (Eq 65):** `H_lA = L2 − LEVCN·CN`, `LEVCN = gp3_lev_cn_ = shells[1]`. Sign
  verified; CN = `es2_energy.coordination`, **verified correct** (matches the binary's setgab CN,
  HF 0.307 vs 0.3062 extracted).
- **Offdiagonals confirmed right:** the new `mfx.exchange_fock` (variational F_x, Euler-verified) +
  `fock_reconstruct.py` isolate the binary's H0 — its offdiagonals are within ~0.01 of ours, so the
  poly/kdiat *structure* is right.

THE REMAINING PIECE — the polar overshoot (HF/H2O), a **compensating multi-error tangle**, not one
wrong parameter:
- Restricting the CN-level to s-shells flips HF −87→+107 — so the p-shell level is a big lever — but
  the CN and LEVCN are the binary's *exact* values, so it's not a wrong CN.
- The basis EEQBC charge terms (SI Eq 28, `gp3_kqvszp_ = l8[2,3,4]`) were tested (`GXTB_BASISQ=1`)
  and **ruled out** — HF barely moves, H2O/NH3 worsen.
- Read: v1's bare baseline is already off for polar (HF −21 before any v2), and the *correct*
  CN-level+Π then over-amplify it. The ledger can't disentangle this (every knob moves several
  molecules at once).

**The clean way forward, and its blocker:** read the binary's *bare* H0 per element and subtract.
`fock_reconstruct.py` gives the exact converged Fock; `H0 = F_ao − ACP − F_x − F_ES`. But the
binary's exchange Fock is **non-variational** (a shortcut potential; `scf_h2.fock_x`), while our
`mfx.exchange_fock` is the exact variational one — so the isolation is clean off-diagonally but the
diagonal (the levels, which is exactly what the polar tangle needs) carries that contamination
(H2 isolates to +0.107). **So the true next step is decoding the binary's non-variational exchange
potential** — then `F_ao − ACP − F_x_nonvar − F_ES = H0` exactly and the tangle reads off directly.
Alternatively, find a pre-SCF bare-Hamiltonian diagonalize (none found so far in the 3 DSYEVR calls).
Also still open and independent: the p-poly *source* for CH4/NH3 (per-element `gp3_pln` overshoots
them; a smaller global helps but that's tuning — resolve via the same clean H0 readout).

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
