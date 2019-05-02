#ifndef ERROR_BOUND_ROUTINES_H
#define ERROR_BOUND_ROUTINES_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include "armadillo"
#include "propagators.hpp"
#include "constants.hpp"

using namespace std;
using namespace arma;
using namespace mgritestimate;

void get_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<double> **lambda, Col<double> *&estimate);

void get_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<cx_double> **lambda, Col<double> *&estimate);

void get_error_l2_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<double> **lambda, Col<double> *&estimate);

void get_error_l2_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<cx_double> **lambda, Col<double> *&estimate);

void get_residual_l2_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<double> **lambda, Col<double> *&estimate);

void get_residual_l2_propagator_bound(const int bound, int theoryLevel, const int cycle, const int relax, Col<int> numberOfTimeSteps, Col<int> coarseningFactors, Col<cx_double> **lambda, Col<double> *&estimate);

double get_error_l2_tight_twogrid_upper_bound(int r, Col<double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_tight_twogrid_upper_bound(int r, Col<cx_double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_tight_twogrid_lower_bound(int r, Col<double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_tight_twogrid_lower_bound(int r, Col<cx_double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_sqrt_expression_upper_bound(int cycle, int r, Col<double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_sqrt_expression_upper_bound(int cycle, int r, Col<cx_double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_sqrt_expression_approximate_rate(int r, Col<double> lambda, Col<int> N, Col<int> m, int theoryLevel);

double get_error_l2_sqrt_expression_approximate_rate(int r, Col<cx_double> lambda, Col<int> N, Col<int> m, int theoryLevel);

void get_samples_index_range(int numberOfSamples, int &samplesRankStartIdx, int &samplesRankStopIdx);

void communicateBounds(Col<double> *&estimate, int numberOfSamples, int samplesRankStartIdx, int samplesRankStopIdx);

#endif