# Paper-vs-binary divergences — the registry

Where the SI (ChemRxiv preprint, our canonical extractions + Richard's screenshots) and the
v1.1 binary (the oracle this build reproduces) demonstrably disagree. Each entry: the evidence,
a MATHEMATICAL analysis of which side is "correct" in what sense, and the implementation
verdict. House rule: the BINARY is the authority for our build (its parameters were fitted
against itself), but every divergence is recorded so (a) upstream can be notified/errata
checked, (b) our Fortran port can offer the theory-consistent variant as an OPTION where one
exists.

---

## D1 — Eq. 55: the repulsion exponent combination rule

- **Paper** (both extractions agree): α_AB = α_A·α_B/(α_A + α_B) — the Gaussian-product
  "reduced exponent" (homonuclear: α/2).
- **Binary** (measured, repulsion gate 4.75e-9 over 13 cases; the halved form refuted 150× by
  the H2 tail): α_AB = 2·α_A·α_B/(α_A + α_B) — the harmonic mean (homonuclear: α).
- **Mathematical analysis**: the reduced exponent is the textbook object for two INTERACTING
  GAUSSIAN charge distributions (it appears in erf-kernel electrostatics), so the paper's form
  has the cleaner pedigree — but the repulsion kernel here is exp(−α_AB·(R+R0)^1.5), and with
  exponent 1.5 (deliberately between Slater and Gaussian, per the SI) NEITHER combination is
  derivable from first principles: it is a parameterization choice, and the fitted α0 values
  only have meaning WITH the rule they were fitted under. Both means are monotone, symmetric,
  and homogeneous degree 1 — no self-consistency test can separate them; only the data can.
- **Verdict**: binary correct-by-construction; paper most plausibly dropped a factor 2 in
  typesetting. No physical error on either side — a definitional mismatch. OUR implementation:
  harmonic (gated); the constant lives in derived-constants.json.

## D2 — the onsite third-order kernel: harmonic vs degree-2

- **Paper**: Eq. 132's onsite branch = 1/(½(1/U + 1/U')) (harmonic, degree 1 in U); the whole
  third order is FRAMED as the charge-derivative of the second-order kernel (Eq. 126).
- **Binary** (measured: the equal-U u² law exact for C/O across ±1/±2; the 2-D grids of
  campaign v2): onsite E3 = c1·q_A·Σ q_l q_l' (U_l+U_l')²/4 — DEGREE 2 in U.
- **Mathematical analysis**: strict DFTB3-style consistency demands E(3) = 1/6 Σ qqq ∂γ⁽²⁾/∂q.
  With γ⁽²⁾_onsite = harmonic(U(q)) and Γ ≡ ∂U/∂q, the equal-U onsite derivative is
  ∂γ/∂U|equal = ½ — degree ZERO in U — times Γ. So: (a) if Γ were constant, theory says
  U-independent onsite E3 — matches NEITHER paper (degree 1) nor binary (degree 2);
  (b) the binary's degree-2 equals the OFFSITE τ formula's leading (U+U')²-structure (Eq. 133a)
  with its R·exp(−kR) factor collapsed to a per-element constant — i.e., the binary looks like
  it applies the offsite machinery at R→0 with a regularized distance factor. NOTE the paper's
  own framing already abandons strict derivative-consistency by introducing a separate γ⁽³⁾.
- **Verdict**: neither side is theory-pure; the binary's form is what the parameters mean.
  OUR implementation: the measured degree-2 form (c1 table). The Fortran port MAY offer a
  "derivative-consistent" theory mode as an option, clearly labelled non-reproducing.

## D3 — Eq. 84: the chemical-potential CN-dependence

- **Paper** (screenshot-canonical): µ = µ⁰(1 + k⁽¹⁾CN·CN_A) — plain CN.
- **Binary** (measured at machine precision via restart-density FD): the coupling tracks
  √CN(mean-rc) at ratio 1.07–1.21 across R = 2..4 while plain CN in either rc convention is off
  by FACTORS. Structure: µ = µ⁰(1 + k·√CN)-family.
- **Mathematical analysis**: CN-scaling of µ is empirical either way (no derivation exists),
  but FAMILY CONSISTENCY strongly favors √CN: the method uses √CN in the repulsion's α(CN)
  (Eq. 56) and in the basis q_eff (Eq. 28's b·√CN term). A √ dropped in typesetting is the
  economical explanation. OPEN: the measured amplitude is ~1.49× the file's k1 slot (suspicious
  near-3/2) with a smooth 6% shape drift — under separation via multi-partner stretches (HF,
  HBr: different pair radii, same H constants).
- **Verdict**: binary (√CN), amplitude question open — entry to be updated by the separation
  round.

---

*Every verdict is reproducible: D1 by `prototype/repulsion.py` (gate), D2 by
`prototype/chart_onsite2.py` + the ion scans, D3 by the restart-precision FD scripts (raw data
in the plan's status log). If upstream publishes source or errata, diff against this registry
first.*

## D4 — Eq. 149: the MFX exchange kernel's Hubbard placement

- **Paper**: γ^MFX = [α + (1−α)erf(ωR)] / (R + favg(U^MFX)·exp(−(k1+k2R)·R)) — the shell
  Hubbards (U^MFX = the L5 row, the only shell-wise row left) in the DENOMINATOR screening:
  |Ex| would shrink as L5 grows.
- **Binary** (measured): ONSITE Ex is EXACTLY linear-increasing in L5 (dEx/dL5 = −0.049268
  constant to six digits over a 16× range — and equal to the master function c_x = s/9.59);
  OFFSITE Ex saturates in L5 (F2 pσ at R=6: L5→0 limit +0.043 of +0.077 pristine, linear-fit
  residual 9.4e-3) — L5 in the numerator AND denominator.
- **Mathematical analysis**: a denominator-only U cannot produce onsite linearity through
  zero intercept (measured exact) under any parameter choice; the onsite channel is the
  separate Sec-1.16 correction with a hardcoded scale (no global moves the atom's Ex — all
  20 give exact zeros), and the offsite kernel carries U multiplicatively with a U-dependent
  screening. The paper's Eq. 149 most plausibly describes an earlier form, or the U-placement
  moved between fit generations; the printed form cannot reproduce the binary.
- **Verdict**: binary (numerator-U onsite-linear + saturating offsite); the measured 17-point
  γ_off(R) curve is the implementation target. CONSEQUENCE: the L5-Euler X-subtraction used
  by the off-diagonal decomposition is exact ONLY onsite/for the tangent part — its offsite
  saturation remainder is what the campaign called "the pzpz object", now reinterpreted.

## D4 addendum — the THEORY adjudication (per the standing request: compare the math to physics)

- The Pariser/Klopman–Ohno tradition defines the exchange-type kernel's onsite limit as the
  chemical hardness: γ(0) ≈ U — the kernel MUST grow with U. The binary's measured onsite
  exchange (exactly linear-increasing in L5, slope = the master function c_x) obeys this;
  the paper's printed denominator-U form (γ(0) = α/U, shrinking with hardness) is
  physically inverted. THEORY VERDICT: the binary implements the physically correct
  convention; Eq. 149 as printed is a misprint, not an alternative model.

## T1 — a binary-vs-THEORY classification (not paper-vs-binary): the non-variational E/F pairing

- Measured: the binary's exchange ENERGY is the exact Mulliken 4-index form (validated on
  four density manifolds to all digits), while its Fock is provably NOT that energy's
  gradient (machine-exact numerical gradients fail the measured Fock on every manifold; the
  population-potential skeleton fits instead). This is the known tight-binding "Mulliken
  shift" shortcut: E from one functional, F from a population approximation.
- CONSEQUENCE (theory-named): the SCF converges to a stationary point of neither
  functional; analytic gradients become inconsistent unless derived against the
  implemented F. The SI's every-section "nuclear gradients will be part of a future
  version" is CONSISTENT with the authors knowing this. Our replication reproduces the
  pairing bug-compatibly (E-formula and F-formula separately), and our future gradient
  work must differentiate the IMPLEMENTED Fock, not the energy.
- Eliminations from the theory tests (fifty-ninth push): the measured short pieces are NOT
  the (true-gradient minus shortcut) difference (ratio drifts 170x across the stretch), so
  the binary's F is not any fixed blend of the two either.
