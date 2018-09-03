debug ?= no
ifeq ($(debug), yes)
    CXXFLAGS = -DDEBUG -g
else
    CXXFLAGS = -DNDEBUG -O2
endif 

default:
	mpic++ src/error_bound_routines.cpp src/sampling_routines.cpp src/rk_routines.cpp src/io_routines.cpp src/propagators.cpp src/operators.cpp main.cpp -o main $(CXXFLAGS) -std=c++11 -Isrc/ -larmadillo

doxygen:
	doxygen doc/doxygen/Doxyfile

load-packages:
	make load-armadillo
	make load-matplotlib

unload-packages:
	make unload-armadillo
	make unload-matplotlib

load-armadillo:
	spack load armadillo

unload-armadillo:
	spack unload armadillo

load-matplotlib:
	spack load py-matplotlib

unload-matplotlib:
	spack unload py-matplotlib