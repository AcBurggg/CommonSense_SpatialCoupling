#!/bin/bash
# Runs the experiment grid used in REPORT.md. Each line of results/*.jsonl is one configuration.
# Usage: bash run_batch.sh <batch-name>
set -u
cd "$(dirname "$0")"
OUT=results
mkdir -p "$OUT"
SIM=./sc_sim
P="n=1000000 d=10000 k=8"   # paper scale: |B| = 1.01M, |B \ A| = 10K, k = 8 (Fig. 2a, first group)
BP="bp_iters=1000 bp_rounds=1"

run() { echo "+ $*" >&2; "$SIM" "$@" >> "$OUT/$BATCH.jsonl"; }

BATCH=${1:-paper}
case "$BATCH" in
paper_uniform)
  for r in 7.6 7.8 8.0; do run $P ratio=$r matrix=uniform decoder=mp trials=16; done
  for r in 7.2 7.4 7.6 7.8; do run $P ratio=$r matrix=uniform decoder=peel trials=16; done
  for r in 5.6 5.8 6.0 6.2; do run $P ratio=$r matrix=uniform decoder=bp trials=16 $BP; done
  ;;
paper_sc)
  for wf in 0.02 0.1; do
    for r in 4.6 4.8 5.0; do run $P ratio=$r matrix=sc wfrac=$wf offsets=tle decoder=peel trials=16; done
    for r in 4.2 4.4 4.6; do run $P ratio=$r matrix=sc wfrac=$wf offsets=tle decoder=bp trials=8 $BP; done
  done
  for r in 4.2 4.4 4.6; do run $P ratio=$r matrix=sc wfrac=0.05 offsets=uni decoder=bp trials=8 $BP; done
  for r in 7.8 8.4; do run $P ratio=$r matrix=sc wfrac=0.05 offsets=tle decoder=mp trials=8; done
  ;;
k_sweep)  # smaller scale (|B| = 101K, |B \ A| = 1K) to locate the best k for SC + BP
  for k in 6 10 12; do
    case $k in 6) RS="5.0 5.5 6.0";; 10) RS="4.0 4.5 5.0";; 12) RS="3.5 4.0 4.5";; esac
    for r in $RS; do run n=100000 d=1000 k=$k ratio=$r matrix=sc wfrac=0.1 offsets=tle decoder=bp trials=16 $BP; done
    case $k in 6) RS="7.5 8.0 8.5";; 10) RS="7.5 8.0 8.5";; 12) RS="8.0 8.5 9.0";; esac
    for r in $RS; do run n=100000 d=1000 k=$k ratio=$r matrix=uniform decoder=mp trials=16; done
  done
  ;;
nd21)  # second operating point: |B \ A| / |B| ~ 1/21 (paper: k = 8, l = 5.1 d)
  Q="n=100000 d=5000 k=8"
  for r in 4.8 5.1 5.4; do run $Q ratio=$r matrix=uniform decoder=mp trials=16; done
  for r in 4.2 4.5 4.8; do run $Q ratio=$r matrix=uniform decoder=bp trials=16 $BP; done
  for r in 3.0 3.3 3.6; do run $Q ratio=$r matrix=sc wfrac=0.1 offsets=tle decoder=bp trials=16 $BP; done
  ;;
esac

# Cheap: exact peeling only (no BP), paper scale, k and window sweeps.
case "$BATCH" in
peel_scan)
  for k in 4 6 8; do
    case $k in 4) RU="6.8 7.0 7.2 7.4"; RS="5.6 5.8 6.0 6.2 6.4";; 6) RU="6.6 6.8 7.0 7.2"; RS="4.6 4.8 5.0 5.2 5.4";; 8) RU="7.0 7.2 7.4 7.6"; RS="4.4 4.6 4.8 5.0 5.2";; esac
    for r in $RU; do run n=1000000 d=10000 k=$k ratio=$r matrix=uniform decoder=peel trials=16 threads=2; done
    for wf in 0.02 0.05 0.1; do
      for r in $RS; do run n=1000000 d=10000 k=$k ratio=$r matrix=sc wfrac=$wf offsets=tle decoder=peel trials=16 threads=2; done
    done
  done
  ;;
esac
