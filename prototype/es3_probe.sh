#!/bin/bash
# es3_probe.sh -- synthetic-input injection probe of set3espot_ (0x53aef0), round 1 (HF).
# Richard's technique: break at the CONVERGED call (hit 9 of 9 for HF), OVERWRITE one input
# channel with synthetic values (zeros / units / doublings), read the output array.
#
# ROUND-1 FINDINGS (2026-07-18, all reproducible by re-running this script):
#  ARGS   rdx -> the LIVE shell-charge array [q_Hs,q_Fs,q_Fp]=[+0.41205,-0.09335,-0.31869]
#         r8  -> the LIVE atomic-charge array [+0.41205,-0.41205]
#         rsi -> per-shell table of per-ELEMENT float32 constants: H -0.35, F -0.22
#                (H2O: O -0.6) -- compiled-in, NOT in the ASCII param file, NOT k_l*row0[6]
#         rdi -> [-0.17645,-0.07047,-0.75003] unidentified (partial-sum suspect)
#         rcx -> OUTPUT V3(3)
#  V3 STRUCTURE (clean, read at the routine boundary):
#         V3(q=0)=0 exactly; V3(2q)=4*V3 BIT-EXACT -> pure homogeneous QUADRATIC in the
#         shell charges; BIT-INDEPENDENT of the r8 atomic-charge argument (the formula's
#         q_A must be built internally as sum_l q_lA).
#  TENSOR (unit+pair injections; FULLY permutation-symmetric -- 9 exact coincidences --
#         and reproduces the converged V3 to 7 digits):
#         T111=0.0617588 T112=0.0082218 T113=0.0128798 T122=0.0051680 T123=0.0066319
#         T133=0.0080959 T222=0.1650071 T223=0.1019432 T233=0.0598809 T333=0.0388203
#         (V3_l = sum_ij T_lij q_i q_j; shells 1=H_s 2=F_s 3=F_p)
#  ENERGY SIDE (open): Euler E3=(1/3)q.V3 = -0.0001815 Eh vs printed-remainder +0.0000243
#         -> the printed ES2+3 is NOT E2_solid + this cubic: the Fock V3 is plausibly
#         NON-VARIATIONAL wrt the energy (precedent: the binary's exchange Fock).
#  WALL (diagnosed): bumping the V3 OUTPUT moves EVERY P-dependent printed term (ES1/ES2+3/
#         multipole/Ex/electronic; repulsion frozen) -> a post-convergence Fock/density
#         rebuild consumes V3, so output-side bumps read dE/dV3 THROUGH density relaxation
#         (the pass-52 wall). Measured relaxed weights dES2+3/dV3 = [-0.3308,+0.0630,+0.2677]
#         (+-0.01 central differences; curvature ~0.4/unit^2 explains the +1.0-bump numbers).
#  NEXT:  match T against SI Eq 131-133 candidates with the known gamma2/U2 diag
#         [0.59406,1.22477,0.59406] + the rsi table; same battery on H2O (4 shells) for
#         cross-molecule validation; identify rdi; source the rsi float32 table's DAT_
#         address from the set3espot_ body in the Ghidra decompile.
EXE=${GXTB_V1_EXE:-/opt/gxtb-v1/gxtb}

run_case() {
  local name="$1" post="$2"; shift 2
  local D
  D=$(mktemp -d)
  cd "$D" || exit 1
  printf '$coord\n 0.0 0.0 0.0 h\n 0.0 0.0 1.733 f\n$end\n' > coord
  {
    echo 'break *0x53aef0'
    echo 'ignore 1 8'
    echo 'run'
    for cmd in "$@"; do echo "$cmd"; done
    echo 'set $out = $rcx'
    echo 'finish'
    if [ -n "$post" ]; then echo "$post"; fi
    echo 'echo \n=== V3POST ===\n'
    echo 'x/3fg $out'
    echo 'continue'
  } > cmds
  echo "##### CASE $name"
  gdb -batch -x cmds --args "$EXE" -c coord 2>/dev/null \
    | grep -E "=== V3POST ===|^0x7f|ES1 |ES2|ES multipole|ES total|Ex \(|electronic|nuclear repulsion"
  cd / && rm -rf "$D"
}

run_case baseline ''
run_case rdx_zero   '' 'set {double}($rdx)=0' 'set {double}($rdx+8)=0' 'set {double}($rdx+16)=0'
run_case r8_zero    '' 'set {double}($r8)=0' 'set {double}($r8+8)=0'
run_case rdx_double '' \
  'set {double}($rdx)=0.8240905807046559' \
  'set {double}($rdx+8)=-0.1867070600559595' \
  'set {double}($rdx+16)=-0.6373839974858573'
run_case B_e1  '' 'set {double}($rdx)=1' 'set {double}($rdx+8)=0' 'set {double}($rdx+16)=0'
run_case B_e2  '' 'set {double}($rdx)=0' 'set {double}($rdx+8)=1' 'set {double}($rdx+16)=0'
run_case B_e3  '' 'set {double}($rdx)=0' 'set {double}($rdx+8)=0' 'set {double}($rdx+16)=1'
run_case B_e12 '' 'set {double}($rdx)=1' 'set {double}($rdx+8)=1' 'set {double}($rdx+16)=0'
run_case B_e13 '' 'set {double}($rdx)=1' 'set {double}($rdx+8)=0' 'set {double}($rdx+16)=1'
run_case B_e23 '' 'set {double}($rdx)=0' 'set {double}($rdx+8)=1' 'set {double}($rdx+16)=1'
