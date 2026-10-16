#!/bin/bash
# Bidirectional experiments (bidir_sim). Usage: bash run_bidir.sh <batch>
set -u
cd "$(dirname "$0")"
mkdir -p results
B=${1:-small}
run() { echo "+ $*" >&2; ./bidir_sim "$@" >> "results/bidir_$B.jsonl"; }
case "$B" in
small)  # |A cap B| = 100K, |A \ B| = |B \ A| = 1K, k = 5 (paper's bidirectional k); ratio = l / |A delta B|
  Q="n=100000 dA=1000 dB=1000 k=5 trials=8"
  for r in 7.5 8.0 8.5 9.0; do run $Q ratio=$r matrix=uniform decoder=mp; done           # paper protocol, 10 rounds max
  for r in 5.0 5.5 6.0 7.0; do run $Q ratio=$r matrix=uniform decoder=bp; done
  for r in 5.0 5.5 6.0; do run $Q ratio=$r matrix=uniform decoder=bp max_rounds=40; done
  for wf in 0.05 0.1 0.25; do for r in 4.0 5.0 6.0; do run $Q ratio=$r matrix=sc wfrac=$wf offsets=tle decoder=bp max_rounds=40; done; done
  for r in 8.5 10.0; do run $Q ratio=$r matrix=sc wfrac=0.1 offsets=tle decoder=mp; done
  ;;
esac

case "$B" in
small_v3)  # same seeds as `small`, now also reporting best-of coders incl. residual diffs
  Q="n=100000 dA=1000 dB=1000 k=5 trials=8"
  for r in 8.5 9.0; do run $Q ratio=$r matrix=uniform decoder=mp; done
  for r in 6.0 7.0; do run $Q ratio=$r matrix=uniform decoder=bp; done
  for wf in 0.05 0.1; do for r in 5.5 6.0 7.0; do run $Q ratio=$r matrix=sc wfrac=$wf offsets=tle decoder=bp max_rounds=40; done; done
  ;;
esac

case "$B" in
asym)  # asymmetric differences: |A \ B| = 200, |B \ A| = 2000 (Bob decodes first, with light noise)
  Q="n=100000 dA=200 dB=2000 k=5 trials=8"
  for r in 4.0 5.0 6.0 7.0; do run $Q ratio=$r matrix=uniform decoder=mp; done
  for r in 3.0 4.0 5.0; do run $Q ratio=$r matrix=uniform decoder=bp max_rounds=40; done
  for r in 2.5 3.0 3.5 4.0; do run $Q ratio=$r matrix=sc wfrac=0.05 offsets=tle decoder=bp max_rounds=40; done
  ;;
esac

case "$B" in
asym2)  # fill-in points for `asym`
  Q="n=100000 dA=200 dB=2000 k=5 trials=8"
  for r in 8.0 9.0 10.0 11.0; do run $Q ratio=$r matrix=uniform decoder=mp; done
  run $Q ratio=4.5 matrix=uniform decoder=bp max_rounds=40
  for r in 4.5 5.0; do run $Q ratio=$r matrix=sc wfrac=0.05 offsets=tle decoder=bp max_rounds=40; done
  ;;
esac
