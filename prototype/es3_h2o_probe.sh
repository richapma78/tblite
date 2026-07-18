#!/bin/bash
# H2O injection battery at set3espot_ converged call (hit 9/9). 4 shells (O_s,O_p,H_s,H_s).
EXE=/opt/gxtb-v1/gxtb

mkcoord_h2o() {
  python3 -c "
import math
B = 1.8897261254578281
ang = 104.5*math.pi/180; r = 0.9572*B
s, c = math.sin(ang/2), math.cos(ang/2)
print('\$coord')
print(' 0.0 0.0 0.0 o')
print(f' {r*s:.10f} 0.0 {r*c:.10f} h')
print(f' {-r*s:.10f} 0.0 {r*c:.10f} h')
print('\$end')
" > coord
}

run_case() {
  local name="$1"; shift
  local D
  D=$(mktemp -d)
  cd "$D" || exit 1
  mkcoord_h2o
  {
    echo 'break *0x53aef0'
    echo 'ignore 1 8'
    echo 'run'
    for cmd in "$@"; do echo "$cmd"; done
    echo 'set $out = $rcx'
    echo 'finish'
    echo 'echo \n=== V3POST ===\n'
    echo 'x/4fg $out'
    echo 'continue'
  } > cmds
  echo "##### CASE $name"
  gdb -batch -x cmds --args "$EXE" -c coord 2>/dev/null \
    | grep -E "=== |^0x7f|ES2"
  cd / && rm -rf "$D"
}

# census: all argument arrays at the converged call
run_case census \
  'echo === rdi4 ===\n' 'x/4fg $rdi' \
  'echo === rsi4 ===\n' 'x/4fg $rsi' \
  'echo === rdx4 ===\n' 'x/4fg $rdx' \
  'echo === rcxpre4 ===\n' 'x/4fg $rcx' \
  'echo === r8_3 ===\n' 'x/3fg $r8'

# structure tests
run_case rdx_zero \
  'set {double}($rdx)=0' 'set {double}($rdx+8)=0' 'set {double}($rdx+16)=0' 'set {double}($rdx+24)=0'
run_case rdx_double \
  'set var *(double*)($rdx)   = 2 * *(double*)($rdx)' \
  'set var *(double*)($rdx+8) = 2 * *(double*)($rdx+8)' \
  'set var *(double*)($rdx+16)= 2 * *(double*)($rdx+16)' \
  'set var *(double*)($rdx+24)= 2 * *(double*)($rdx+24)'
run_case r8_zero \
  'set {double}($r8)=0' 'set {double}($r8+8)=0' 'set {double}($r8+16)=0'

# units
for i in 0 1 2 3; do
  cmds=()
  for j in 0 1 2 3; do
    v=0; [ "$i" = "$j" ] && v=1
    cmds+=("set {double}(\$rdx+$((j*8)))=$v")
  done
  run_case "e$((i+1))" "${cmds[@]}"
done

# pairs
for p in "0 1" "0 2" "0 3" "1 2" "1 3" "2 3"; do
  set -- $p
  a=$1; b=$2
  cmds=()
  for j in 0 1 2 3; do
    v=0; { [ "$j" = "$a" ] || [ "$j" = "$b" ]; } && v=1
    cmds+=("set {double}(\$rdx+$((j*8)))=$v")
  done
  run_case "e$((a+1))$((b+1))" "${cmds[@]}"
done
