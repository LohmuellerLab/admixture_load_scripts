#! /bin/bash

for model in "model0" "model4"; do

(
  echo "rep,r,splittime,generation,p1Fraction,p2Fraction,FST_bhatia,FST_hudson,FST_mean"
  
  for s in 100 250 500 1000 2500 5000 10000 15000 20000 25000 30000 35000 40000; do
  
  for r in "1e-9"; do
    for i in {1..200}; do
      cat ./sims/${model}/split_${s}/slim_h0.0_r${r}_${i}.csv | sed -e '1,/#generation/ d' |  awk -v var="$i,$r,$s," '{print var $0}'
    done
  done
  
  done
) > "${model}_split".csv

done