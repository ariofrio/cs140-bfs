#!/bin/sh
if [ -z "$PROCS" ]; then
  PROCS="32 16 8 4 2 1"
fi

if [ -z "$ALGOS" ]; then
  ALGOS="1 2 3"
fi

for p in $PROCS; do
  for a in $ALGOS; do
    echo "Running algorithm $a with $p workers"
    time ../parallel$a $1 $2 -v -cilk_set_worker_count=$p \
      2>&1 | tee results/$(basename $1).from$2.workers$p.$a
  done
done

