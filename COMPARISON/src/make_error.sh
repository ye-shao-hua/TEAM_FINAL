#!/bin/bash
./bin/error.out ./result/real_cap/input_diag.txt ./result/calculated_cap/input_diag.txt ./result/error/input_diag.txt 48
./bin/error.out ./result/real_cap/input_plate2l.txt ./result/calculated_cap/input_plate2l.txt ./result/error/input_plate2l.txt 32
./bin/error.out ./result/real_cap/input_plate3l.txt ./result/calculated_cap/input_plate3l.txt ./result/error/input_plate3l.txt 32
./bin/error.out ./result/real_cap/input_stack.txt ./result/calculated_cap/input_stack.txt ./result/error/input_stack.txt 48

for i in `seq 1 3`
do
  ./bin/error.out ./result/real_cap/PLATE2L_metal${i}.txt ./result/calculated_cap/PLATE2L_metal${i}.txt ./result/error/PLATE2L_metal${i}.txt 324
done

for i in `seq 1 2`
do 
  for j in `seq $((i+1)) 3`
  do 
    ./bin/error.out ./result/real_cap/PLATE3L_metal${i}_metal${j}.txt ./result/calculated_cap/PLATE3L_metal${i}-metal${j}.txt ./result/error/PLATE3L_metal${i}-metal${j}.txt 324
    ./bin/error.out ./result/real_cap/STACK3L_metal${i}_metal${j}.txt ./result/calculated_cap/STACK3L_metal${i}-metal${j}.txt ./result/error/STACK3L_metal${i}-metal${j}.txt 144
  done
done
