#!/bin/bash
for i in `seq 1 48`
do
	./../bin/solvor -i ./../Cases/input_diag/BEM_INPUT_${i}_* -o ./../result/input_diag/out_${i}.txt
	./../bin/solvor -i ./../Cases/input_stack/BEM_INPUT_${i}_* -o ./../result/input_stack/out_${i}.txt
done

for i in `seq 1 32`
do
	./../bin/solvor -i ./../Cases/input_plate2l/BEM_INPUT_${i}_* -o ./../result/input_plate2l/out_${i}.txt
	./../bin/solvor -i ./../Cases/input_plate3l/BEM_INPUT_${i}_* -o ./../result/input_plate3l/out_${i}.txt
done

for j in `seq 1 3`
do
       for k in `seq 1 324`
       do
	./../bin/solvor -i ./../Cases/PLATE2L/SUB-metal${j}_PLATE2L/input/BEM_INPUT_${k}_* -o ./../result/PLATE2L/m${j}/out_${l}.txt

        done
done

for k in `seq 1 324`
do
	#echo "$(expr 2 \* $k)"
        ./../bin/solvor -i ./../Cases/PLATE3L/SUB-metal1-metal2_PLATE3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/PLATE3L/m1-m2/out_${k}.txt
        ./../bin/solvor -i ./../Cases/PLATE3L/SUB-metal1-metal3_PLATE3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/PLATE3L/m1-m3/out_${k}.txt
        ./../bin/solvor -i ./../Cases/PLATE3L/SUB-metal2-metal3_PLATE3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/PLATE3L/m2-m3/out_${k}.txt
done

for k in `seq 1 144`
do
	#echo "$(expr 2 \* $k)"
        ./../bin/solvor -i ./../Cases/STACK3L/SUB-metal1-metal2_STACK3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/STACK3L/m1-m2/out_${k}.txt
        ./../bin/solvor -i ./../Cases/STACK3L/SUB-metal1-metal3_STACK3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/STACK3L/m1-m3/out_${k}.txt
        ./../bin/solvor -i ./../Cases/STACK3L/SUB-metal2-metal3_STACK3L/input/BEM_INPUT_$(expr 2 \* $k)_* -o ./../result/STACK3L/m2-m3/out_${k}.txt
done

