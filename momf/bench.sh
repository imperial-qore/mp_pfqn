#!/bin/bash
# Wall-clock comparison of the exact rational solver against the
# multi-modular one, as the population grows.
BASE=$1; shift
SCALES="$@"
R=$(dirname $0)/..
printf "%6s %10s %12s %12s %12s %8s\n" "Ntot" "digits(G)" "mom(s)" "mommod-j1" "mommod-j8" "speedup"
for s in $SCALES; do
  f=$(mktemp /tmp/bench_XXXX.qn)
  awk -v s=$s 'NR==2{for(i=1;i<=NF;i++)$i=$i*s} {print}' $BASE > $f
  t0=$(date +%s.%N); $R/bin/mom    $f -e > /tmp/bench_g.txt 2>/dev/null; t1=$(date +%s.%N)
  dg=$(head -1 /tmp/bench_g.txt | wc -c)
  t2=$(date +%s.%N); $R/bin/mommod $f -e -j 1 > /dev/null 2>&1;         t3=$(date +%s.%N)
  t4=$(date +%s.%N); $R/bin/mommod $f -e -j 8 > /dev/null 2>&1;         t5=$(date +%s.%N)
  a=$(echo "$t1-$t0"|bc); b=$(echo "$t3-$t2"|bc); c=$(echo "$t5-$t4"|bc)
  sp=$(echo "scale=2; $a/$c"|bc)
  n=$(awk 'NR==2{t=0;for(i=1;i<=NF;i++)t+=$i;print t}' $f)
  printf "%6s %10s %12.3f %12.3f %12.3f %8s\n" $n $dg $a $b $c $sp
  rm -f $f
done
