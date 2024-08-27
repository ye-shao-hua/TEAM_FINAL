CC=g++
BIN=./bin
all: solvor

solvor: src/transfer.cpp ./lib/libcopredict1.a ./lib/libcopredict2.a ./lib/libcopredict3.a ./lib/libcopredict4.a ./lib/libcopredict5.a ./lib/libcopredict6.a ./lib/libcopredict7.a
	${CC} ./src/transfer.cpp -o ./bin/solvor -I include -std=c++23 -I ./etc/tiny-dnn -I ./etc/CLI11/include -lcopredict1 -lcopredict2 -lcopredict3 -lcopredict4 -lcopredict5 -lcopredict6 -lcopredict7 -L./lib -static

./lib/libcopredict1.a: src/copredict1.cpp
	${CC} ./src/copredict1.cpp -c -o ./tmp/copredict1.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict1.a ./tmp/copredict1.o

./lib/libcopredict2.a: src/copredict2.cpp
	${CC} ./src/copredict2.cpp -c -o ./tmp/copredict2.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict2.a ./tmp/copredict2.o

./lib/libcopredict3.a: src/copredict3.cpp
	${CC} ./src/copredict3.cpp -c -o ./tmp/copredict3.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict3.a ./tmp/copredict3.o

./lib/libcopredict4.a: src/copredict4.cpp
	${CC} ./src/copredict4.cpp -c -o ./tmp/copredict4.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict4.a ./tmp/copredict4.o

./lib/libcopredict5.a: src/copredict5.cpp
	${CC} ./src/copredict5.cpp -c -o ./tmp/copredict5.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict5.a ./tmp/copredict5.o

./lib/libcopredict6.a: src/copredict6.cpp
	${CC} ./src/copredict6.cpp -c -o ./tmp/copredict6.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict6.a ./tmp/copredict6.o

./lib/libcopredict7.a: src/copredict7.cpp
	${CC} ./src/copredict7.cpp -c -o ./tmp/copredict7.o -std=c++17 -I ./etc/tiny-dnn -static
	ar rcs ./lib/libcopredict7.a ./tmp/copredict7.o


clean:
	rm -f ${BIN}/*
	rm -f ./lib/*
	rm -f ./tmp/*
