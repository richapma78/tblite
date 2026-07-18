#!/bin/bash
# R-sweep of the ES3 tau kernel + gamma2 (both packed, read at EVERY hit; the last block
# per routine is the converged one). H2 (3 packed) and HF (6 packed).
EXE=/opt/gxtb-v1/gxtb

sweep_point() { # $1=label $2=coord-content $3=npacked
  local D
  D=$(mktemp -d)
  cd "$D" || exit 1
  printf '%s\n' "$2" > coord
  {
    echo 'break *0x53aef0'
    echo 'commands 1'
    echo 'silent'
    echo 'echo TAU\n'
    echo "x/$3fg \$rdi"
    echo 'continue'
    echo 'end'
    echo 'break *0x53ab60'
    echo 'commands 2'
    echo 'silent'
    echo 'echo GAM\n'
    echo "x/$3fg \$rdi"
    echo 'continue'
    echo 'end'
    echo 'run'
  } > cmds
  echo "@@POINT $1"
  gdb -batch -x cmds --args "$EXE" -c coord 2>/dev/null | grep -E "^TAU|^GAM|^0x7f"
  cd / && rm -rf "$D"
}

for R in 1.2 1.4 1.7 2.0 2.4 3.0 3.6 4.5; do
  sweep_point "H2_$R" "\$coord
 0.0 0.0 0.0 h
 0.0 0.0 $R h
\$end" 3
done

for R in 1.4 1.733 2.0 2.5 3.0 4.0; do
  sweep_point "HF_$R" "\$coord
 0.0 0.0 0.0 h
 0.0 0.0 $R f
\$end" 6
done
