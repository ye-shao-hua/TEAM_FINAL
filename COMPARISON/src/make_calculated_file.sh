#!/bin/bash
for i in `seq 1 3`
do 
  n=1
  for j in `seq 1 324`
  do 
    ./bin/readOut.out ./tmp/PLATE2L/m${i}/out_${j}.txt 5 ${n}>> ./result/calculated_cap/PLATE2L_metal${i}.txt
    n=`expr $n + 1`
  done
done

for i in `seq 1 2`
do 
  for j in `seq $((i+1)) 3`
  do
    n=1
    for k in `seq 1 324`
    do
      ./bin/readOut.out ./tmp/PLATE3L/m${i}-m${j}/out_${k}.txt 5 ${n} >> ./result/calculated_cap/PLATE3L_metal${i}-metal${j}.txt
      n=`expr $n + 1`
    done
  done
done

for i in `seq 1 2`
do 
  for j in `seq $((i+1)) 3`
  do
    n=1
    for k in `seq 1 144`
    do
      ./bin/readOut.out ./tmp/STACK3L/m${i}-m${j}/out_${k}.txt 9 ${n} >> ./result/calculated_cap/STACK3L_metal${i}-metal${j}.txt
      n=`expr $n + 1`
    done
  done
done

n=1
for i in `seq 1 48`
do
  ./bin/readOut.out ./tmp/input_diag/out_${i}.txt 9 ${n} >> ./result/calculated_cap/input_diag.txt
  n=`expr $n + 1`
done

n=1
for i in `seq 1 32`
do
  ./bin/readOut.out ./tmp/input_plate2l/out_${i}.txt 6 ${n} >> ./result/calculated_cap/input_plate2l.txt
  n=`expr $n + 1`
done

n=1
for i in `seq 1 32`
do
  ./bin/readOut.out ./tmp/input_plate3l/out_${i}.txt 8 ${n} >> ./result/calculated_cap/input_plate3l.txt
  n=`expr $n + 1`
done

n=1
for i in `seq 1 48`
do
  ./bin/readOut.out ./tmp/input_stack/out_${i}.txt 9 ${n} >> ./result/calculated_cap/input_stack.txt
  n=`expr $n + 1`
done
