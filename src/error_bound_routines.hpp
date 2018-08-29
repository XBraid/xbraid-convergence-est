#ifndef ERROR_BOUND_ROUTINES_H
#define ERROR_BOUND_ROUTINES_H

#include "mpi.h"
#include <iostream>
#include "armadillo"
#include "propagators.hpp"
#include "constants.hpp"

using namespace std;
using namespace arma;
using namespace mgritestimate;

void get_error_propagator_bound(const int bound, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<cx_double> **lambda, Col<double> *&estimate);

#endif