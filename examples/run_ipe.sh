#!/bin/bash
# Minimal example: run standalone IPE for one segment.
# Set the paths below to your environment, then launch with mpiexec/mpirun.
set -e

IPEEXEC=${IPEEXEC:-../install/bin/ipe.x}
FIX_DIR=${FIX_DIR:?directory containing IPE_Grid.h5 and the fix files (wei96*, *.dat, ...)}
IC=${IC:?path to an IPE_State.apex.*.h5 initial condition}
FORCING=${FORCING:?path to wam_input_f107_kp.txt}
RUNDIR=${RUNDIR:-./run}
NPROC=${NPROC:-40}
CDATE=${CDATE:-2023010100}          # must match INITIAL_TIMESTAMP in IPE.inp

mkdir -p "$RUNDIR/output"
cd "$RUNDIR"
ln -sf "$FIX_DIR"/* .
cp "$IC" "IPE_State.apex.${CDATE}00.h5"
cp "$FORCING" wam_input_f107_kp.txt
cp "$(dirname "$0")/IPE.inp.example" IPE.inp     # edit for your run
cp "$IPEEXEC" ./ipe.x

mpiexec -n "$NPROC" ./ipe.x
