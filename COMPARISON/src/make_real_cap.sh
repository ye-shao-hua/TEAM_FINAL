#!/bin/bash
./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal1_PLATE2L/SUB-metal1_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal1.txt
./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal2_PLATE2L/SUB-metal2_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal2.txt
./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal3_PLATE2L/SUB-metal3_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal3.txt
clear
./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal1-metal2_PLATE3L/SUB-metal1-metal2_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal1_metal2.txt
./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal1-metal3_PLATE3L/SUB-metal1-metal3_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal1_metal3.txt
./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal2-metal3_PLATE3L/SUB-metal2-metal3_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal2_metal3.txt
clear
./bin/readStack3l.out ./Cases/STACK3L/SUB-metal1-metal2_STACK3L/SUB-metal1-metal2_STACK3L.tbl.text ./result/real_cap/STACK3L_metal1_metal2.txt
./bin/readStack3l.out ./Cases/STACK3L/SUB-metal1-metal3_STACK3L/SUB-metal1-metal3_STACK3L.tbl.text ./result/real_cap/STACK3L_metal1_metal3.txt
./bin/readStack3l.out ./Cases/STACK3L/SUB-metal2-metal3_STACK3L/SUB-metal2-metal3_STACK3L.tbl.text ./result/real_cap/STACK3L_metal2_metal3.txt
clear
cp ./Cases/input_diag/table_result_diag.txt ./result/real_cap/
cp ./Cases/input_plate2l/table_result_plate2l.txt ./result/real_cap/
cp ./Cases/input_plate3l/table_result_plate3l.txt ./result/real_cap/
cp ./Cases/input_stack/table_result_stack.txt ./result/real_cap/
clear
