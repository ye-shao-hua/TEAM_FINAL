CC=g++
STD=-std=c++11
all:make_result generator execute readOut.cpp make_error

make_error:readOut.cpp execute
	chmod +x ./src/make_error.sh
	./src/make_error.sh
	#clear
readOut.cpp:make_result
	${CC} ./src/readOut.cpp -o ./bin/readOut.out ${STD}
	chmod +x ./src/make_calculated_file.sh
	./src/make_calculated_file.sh
	#clear
execute:generator make_result
	#chmod +x ./src/make_real_cap.sh
	#./src/make_real_cap.sh
	./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal1_PLATE2L/SUB-metal1_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal1.txt
	./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal2_PLATE2L/SUB-metal2_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal2.txt
	./bin/readPlate2l.out ./Cases/PLATE2L/SUB-metal3_PLATE2L/SUB-metal3_PLATE2L.tbl.text ./result/real_cap/PLATE2L_metal3.txt
	#clear
	./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal1-metal2_PLATE3L/SUB-metal1-metal2_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal1_metal2.txt
	./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal1-metal3_PLATE3L/SUB-metal1-metal3_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal1_metal3.txt
	./bin/readPlate3l.out ./Cases/PLATE3L/SUB-metal2-metal3_PLATE3L/SUB-metal2-metal3_PLATE3L.tbl.text ./result/real_cap/PLATE3L_metal2_metal3.txt
	#clear
	./bin/readStack3l.out ./Cases/STACK3L/SUB-metal1-metal2_STACK3L/SUB-metal1-metal2_STACK3L.tbl.text ./result/real_cap/STACK3L_metal1_metal2.txt
	./bin/readStack3l.out ./Cases/STACK3L/SUB-metal1-metal3_STACK3L/SUB-metal1-metal3_STACK3L.tbl.text ./result/real_cap/STACK3L_metal1_metal3.txt
	./bin/readStack3l.out ./Cases/STACK3L/SUB-metal2-metal3_STACK3L/SUB-metal2-metal3_STACK3L.tbl.text ./result/real_cap/STACK3L_metal2_metal3.txt
	#clear
	cp ./Cases/input_diag/table_result_diag.txt ./result/real_cap/input_diag.txt
	cp ./Cases/input_plate2l/table_result_plate2l.txt ./result/real_cap/input_plate2l.txt
	cp ./Cases/input_plate3l/table_result_plate3l.txt ./result/real_cap/input_plate3l.txt
	cp ./Cases/input_stack/table_result_stack.txt ./result/real_cap/input_stack.txt
	#clear
	cp -r ./../NETWORK/result ./tmp
	#clear
make_result:./src/make_result.sh
	make clean
	chmod +x ./src/make_result.sh
	./src/make_result.sh
	#clear
generator:./src/*.cpp
	$(CC)  ./src/readPlate2l.cpp -o ./bin/readPlate2l.out ${STD}
	$(CC)  ./src/readPlate3l.cpp -o ./bin/readPlate3l.out	${STD}
	$(CC)  ./src/readStack3l.cpp -o ./bin/readStack3l.out	${STD}
	${CC}  ./src/error.cpp 			 -o ./bin/error.out 			${STD}
clean:
	rm -f ./lib/*
	rm -rf ./tmp/
	rm -f ./bin/*
	rm -rf ./result
	#clear
