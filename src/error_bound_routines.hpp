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

void get_error_propagator_bound(const int bound, int theoryLevel, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<cx_double> **lambda, Col<double> *&estimate);

double get_sqrt_expression_upper_bound(int r, Col<cx_double> lambda, Col<int> numberOfTimeSteps, Col<int> m, int theoryLevel);

void get_samples_index_range(int numberOfSamples, int &samplesRankStartIdx, int &samplesRankStopIdx);

void communicateBounds(Col<double> *&estimate, int numberOfSamples, int samplesRankStartIdx, int samplesRankStopIdx);

#endif