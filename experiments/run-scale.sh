#!/bin/sh
if [ -z "$PROCS" ]; then
  PROCS="32 16 8 4 2 1"
fi

for p in $PROCS; do
  echo "Running naive linked list with $p workers"
  time ../parallel1 $1 $2 -v -cilk_set_worker_count=$p \
    2>&1 | tee results/$(basename $1).from$2.workers$p.naive
  echo "Running smart linked list with $p workers"
  time ../parallel2 $1 $2 -v -cilk_set_worker_count=$p \
    2>&1 | tee results/$(basename $1).from$2.workers$p.smart
done

