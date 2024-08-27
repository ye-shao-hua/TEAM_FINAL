XX=g++
TARGET=test.out
SRC=test.cpp
STD=-std=c++23

all:test.out run

test.out:test.cpp libdie_manager.a
	${XX} ${SRC} ${STD} -o ${TARGET} -I include -O3 -fopenmp -mavx2 -msse -mfma -ldie_manager -L.
	#${XX} ${SRC} ${STD} -o ${TARGET} -I include -O3 -mavx2 -msse -mfma -ldie_manager -L.

die_manager.o:die_manager.cpp
	${XX} -c die_manager.cpp -I include/ -o die_manager.o -O3

libdie_manager.a:die_manager.o
	ar rcs libdie_manager.a die_manager.o

run: test.out
	./${TARGET}

clean:
	rm ~/ysh/${TARGET}

