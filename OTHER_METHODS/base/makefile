CXXFLAG=-O3 -flto -fopenmp

all:main.out
	clear
	./bin/main.out

%.o:src/%.cpp
	g++ -c -o ./tmp/$@ $< -std=c++11 ${CXXFLAG}

libfem.a:full_matrix.o full_vector.o base_function.o point.o sparsity_matrix.o sparse_matrix.o cell.o fe_q.o hyper_cube.o mapping_q.o triangulation.o
	ar rcs ./lib/$@ ./tmp/*

main.out:main.cpp libfem.a
	g++ $< -lfem -L./lib/ -I./include/ -o ./bin/$@ -std=c++11 ${CXXFLAG} -g

clear:
	rm -rf ./bin/*
	rm -rf ./lib/*
	rm -rf ./tmp/*
