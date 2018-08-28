default:
	mpic++ src/rk_routines.cpp src/io_routines.cpp src/propagators.cpp src/operators.cpp main.cpp -o main -O2 -std=c++11 -Isrc/ -larmadillo
