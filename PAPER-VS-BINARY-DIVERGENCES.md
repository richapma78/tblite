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
