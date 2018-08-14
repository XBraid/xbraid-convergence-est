default:
	mpic++ src/propagators.cpp src/operators.cpp main.cpp -o main -O2 -std=c++11 -Isrc/ -larmadillo
