# The `gp3_` Parameter Symbol Table — the g-xTB binary's Rosetta Stone

**Found 2026-07-18.** The reference g-xTB executable (`/opt/gxtb-v1/gxtb`, an ifort build with
symbols intact) keeps **every** method parameter in a named module array, prefix
`parcom_mp_gp3_*`. The symbol *names* say which column of the ASCII parameter file
(`~/.gxtb`, parsed by `paramfile.py` / `params.py`) feeds which equation of the SI — which
ends the column-guessing the port had been doing by hand.

This is the authoritative, per-element source for the whole method's parameter layer. Where a
prototype value was previously *derived* or *fit*, the matching `gp3_` array is the ground
truth to check it against.

## How to read them

```bash
# list every parameter array
nm /opt/gxtb-v1/gxtb | grep parcom_mp_gp3

# read one live (break AFTER the param file is loaded; setespot is a convenient anchor)
gdb -batch \
  -ex 'break *0x53ab60' -ex run \
  -ex 'x/40fg 0x3be42e0' \        # e.g. gp3_lev_
  --args /opt/gxtb-v1/gxtb -c coord      # any coord file in the cwd
```

- Each array is **per element**, stride **4 doubles** (`s, p, d, f`) unless noted; element 1 =
  H at offset 0, element `Z` at offset `(Z-1) * 32` bytes.
- Values are bit-identical to the ASCII parameter file — the arrays are just the *named,
  column-resolved* view of it.
- The address column below is from this build; re-derive with `nm` if the binary changes.

## The map (SI section in parentheses)

| symbol | addr | = param-file column | meaning |
|---|---|---|---|
| `gp3_lev_` | `0x3be42e0` | `shells[0]` (= `L2`) | H0 base shell level `h_lA` (Eq 64/65). H_s = 0.21619978 |
| `gp3_lev_cn_` | `0x3be5640` | `shells[1]` | H0 shell-level **CN dependence** `k^{H,CN}` (Eq 65: `H_lA = h_lA − k^{H,CN}·CN_A`). H_s = −0.01708. **GE.build drops this** |
| `gp3_poly_` | `0x3bfb5e0` | `l1[7]` (row0[7]) | H0 shell-polynomial element coeff `k^{shp}_A` (Eq 67). H = 0.02236770 |
| `gp3_rad_` | `0x3bfaf60` | `l1[5]` (`rcov_cn`) | the shell-polynomial covalent radius `R^cov` (Eq 67). H = 1.05463349. *Distinct from the CN-count radius* (`es2_energy.RCOV` ≈ 0.665) — the two-radii theme |
| `gp3_h0k_` | `0x3be4fc0` | `l8[0]` | H0 diatomic scaling `kdiat` (GE `kd_sg`). H = 2.90912803 |
| `gp3_pln_` | `0x3c00360` | `shells[2]` | q-vSZP exponent scale `k` (GE `E[z]['k']`). H = 1.09404075 |
| `gp3_glob_` | `0x3be6320` | `g1[0..9] ++ g2[0..9]` | the 20 globals (`kW_s/p = g1[0]/g1[1]`). The 4 `k^{shp,l}` H0 globals are **mixed into these**, not a separate slot — the one H0 piece still unsourced |
| `gp3_gam1_` | — | ES1 Γ | first-order shell Hubbard-derivative table |
| `gp3_gam2_` | `0x3be32c0` | ES2 Γ | second-order table |
| `gp3_gam3_` | `0x3a7eca0` | **ES3 Γ** | third-order table. **CONFIRMED bit-identical to the `rsi` table `set3espot_` loads** (H=−0.35, O=−0.60, F=−0.22) → ES3's `Γ_A` **is** `gp3_gam3_[Z]` |
| `gp3_e1cn_` | `0x3bff660` | ES1 CN-dep | shell-charge CN dependence for ES1 |
| `gp3_g2cn_` | `0x3bfb2a0` | ES2 CN-dep | shell-charge CN dependence for ES2 |
| `gp3_acp_` | `0x3bfd960` | ACP | atomic-correction-potential params |
| `gp3_acpcn_` | `0x3bff320` | ACP CN-dep | " with CN dependence |
| `gp3_rep_` | `0x3bfb920` | repulsion | semi-classical repulsion params |
| `gp3_aesr0_` | `0x3be63c0` | AES `R0` | anisotropic-electrostatics damping table (= `aes.R0`) |
| `gp3_kqvszp_` | `0x3bff9a0` | q-vSZP `k̃` | basis charge-dependence coefficients |
| `gp3_atspin_` | `0x3be22a0` | atomic spin | spin-polarisation params |
| `gp3_at_`, `gp3_fam_`, `gp3_e1cn_` siblings | — | — | further per-element tables (name-decode as needed) |

## Why it matters

- **Sources the whole per-element layer by name** — no more matching a fit against an unnamed
  column. Every `gp3_*` is the ground truth for its term.
- **Unlocked the core Hamiltonian (H0)**: the four per-element pieces of SI Eq 64-67 (`gp3_lev_`,
  `gp3_lev_cn_`, `gp3_poly_`, `gp3_rad_`) are now sourced; only the 4 `k^{shp,l}` globals and the
  σ/π `kdiat` treatment remain. See `data/derived-constants.json` → `h0_electronic_gap_ROADMAP`
  and the banked patch `data/h0_v2_CN_levels_and_poly.patch`.
- **Resolved the ES3 Γ**: `gp3_gam3_` is the `rsi` table `set3espot_` reads, closing a
  standing "unidentified" note.

Full machine-readable map: `data/derived-constants.json` → `gp3_parameter_symbols_ROSETTA`.
Background on the binary and the decompile: the `gxtb-decompilation-asset` memory.
