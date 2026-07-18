#!/bin/bash
# aes_epol_probe.sh -- extract stock xtb_aespot aniso_electro inputs/outputs to decode the
# g-xTB "ES multipole" energy. READS ONLY (injection breaks the buffered-stdout flush;
# see below). Run: bash aes_epol_probe.sh
#
# THE FINDINGS (2026-07-18), all reproducible:
#  * g-xTB calls STOCK xtb_aespot (out-of-line, 0x4cc620=aniso_electro, NOT inlined). The
#    exact energy formula is the stock aniso_electro_cpu (reference xtb-bleed aespot.F90:569):
#      e01 = sum gab3*(charge-dipole)       [gab3 = R^-3 erf-damped kernel]
#      e02 = sum gab5*(charge-quadrupole)   [gab5 = R^-5 kernel]
#      e11 = sum gab5*(dipole-dipole)
#      epol= sum_i [dipKernel(Zi)*|dip_i|^2 + quadKernel(Zi)*|qp_i|^2]   (the SELF term)
#      printed "ES multipole" = e01+e02+e11+epol   (g-xTB folds epol into e)
#  * VERIFIED bit-exact against the binary: gab3, gab5 (0.05260597.., 0.00521734..) and the
#    CAMM dipoles (H2 |dip|=0.28469252) match aes.py EXACTLY. So moments+kernels are right.
#  * aes.py's ENERGY was WRONG in TWO ways: (1) it ADDED spurious g7 (dip-quad) and g9
#    (quad-quad) orders -- STOCK xtb HAS ONLY gab3/gab5, no gab7/gab9; (2) it OMITTED epol.
#    On the homonuclear controls (q=0 kills e01/e02) the whole miss IS epol: H2 -0.698, F2
#    -0.800 mEh = exactly the epol this decode found (H2 epol=+0.672, e11=+1.658, e=+2.329).
#  * epol's kernels are CN-DEPENDENT (the g-xTB modification; stock kernels are constant).
#    Constant-kernel fit over 5 bond lengths floors at ~58 uEh (systematic curve); a
#    CN-linear kernel dipKernel=a+b*CN closes H to 0.7 uEh, F to 2.2 uEh. CN = the internal
#    es2_energy.coordination. (The struct TMultipoleData carries cnShift/valenceCN for this.)
#  * OPEN: the exact per-element base kernels + CN form for H,C,N,O,F. Homonuclear diatomics
#    only give H,F (fitted a,b: H dipK=0.00643-0.00327*CN, quadK=-0.00062-0.00028*CN; F
#    dipK=0.00741+0.01037*CN). C,N,O need another source: the g-xTB SI AES self-energy
#    section, or the aesData struct arrays (dipKernel/quadKernel allocatables), or CAMM
#    self-term injection at HIT 1 done without breaking the flush.
#
# BUFFERED-STDOUT TRAP (Richard's tip, confirmed): aniso_electro is called 3x; the printed
# "ES multipole" is computed at HIT 1 but Fortran block-buffers stdout when piped, so it only
# FLUSHES at normal exit. Zeroing an input at ANY hit corrupts the tail and the buffered line
# is lost -- so array-zeroing injection reads NOTHING. READ (don't inject) at hit 1, or find
# the flush. All values below are READS.
EXE=${GXTB_V1_EXE:-/opt/gxtb-v1/gxtb}

extract() { # $1=sym $2=R -- read moments+energy at hit 1 (the print-feeding call)
  local D; D=$(mktemp -d); cd "$D" || exit 1
  printf '$coord\n 0.0 0.0 0.0 %s\n 0.0 0.0 %s %s\n$end\n' "$1" "$2" "$1" > coord
  cat > cmds <<EOF
break *0x4cc620
run
set \$p5=*(long*)\$r8
set \$p6=*(long*)\$r9
set \$p8=*(long*)(\$rsp+16)
set \$e=*(long*)(\$rsp+40)
printf "DIPM %.12f %.12f\n", *(double*)(\$p5+16), *(double*)(\$p5+40)
printf "QP %.12f %.12f %.12f\n", *(double*)\$p6, *(double*)(\$p6+16), *(double*)(\$p6+40)
printf "GAB5 %.12f\n", *(double*)(*(long*)\$p8+8)
finish
printf "E %.12f\n", *(double*)\$e
EOF
  echo "@@ $1 $2"
  gdb -batch -x cmds --args "$EXE" -c coord 2>/dev/null | grep -E "DIPM|QP |GAB5|^E "
  cd / && rm -rf "$D"
}

for R in 1.2 1.4 1.7 2.0 2.4; do extract h $R; done
for R in 2.0 2.2 2.5 2.668 3.2; do extract f $R; done
