#! /bin/bash

for type in "model0" "model1" "model2" "model3" "model4"; do

    cd burnin/${type}
    for h in 0.0 0.5; do

        cd h_${h}

(
  echo "rep,r,generation,meanFitnessP1,meanFitnessP2,p1Fraction,p1FractionWithin,p1FractionOutside,p2Fraction,p2FractionWithin,p2FractionOutside,p1MeanDel,p1MeanNeu,p1MeanMut,p1IndivDel,p1MeanHomNeu,p1MeanHomDel,p1load_s,p1load_w,p2MeanDel,p2MeanNeu,p2MeanMut,p2IndivDel,p2MeanHomNeu,p2MeanHomDel,p2load_s,p2load_w"
  for r in "1e-6" "1e-7" "1e-8" "1e-9"
  do
    for i in {1..200}
    do
      cat slim_h${h}_r${r}_${i}.csv | sed -n -e '/#generation/,$p' | awk '/#generation,FitnessLarge,FitnessSmall,p1Fraction,p1FractionWithin,p1FractionOutside,p2Fraction,p2FractionWithin,p2FractionOutside,p1MeanDel,p1MeanNeu,p1MeanMut,p1IndivDel,p1MeanHomNeu,p1MeanHomDel,p1load_s,p1load_w,p2MeanDel,p2MeanNeu,p2MeanMut,p2IndivDel,p2MeanHomNeu,p2MeanHomDel,p2load_s,p2load_w/{y=1;next}y' | awk -v var="$i,$r," '{print var $0}'
    done
  done
) > "h_${h}_${type}".csv

        cp "h_${h}_${type}".csv ../../

        cd ..
    
    done

    cd ../..

done