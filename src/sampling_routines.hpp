#ifndef SAMPLING_ROUTINES_H
#define SAMPLING_ROUTINES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"
#include "rk_routines.hpp"

using namespace rkconst;

void sample_complex_plane(const int method,
                          double minreal, double maxreal, double minimag, double maximag,
                          int numberOfRealSamples, int numberOfImagSamples,
                          arma::Col<arma::cx_double> *&z, arma::Col<arma::cx_double> *&lambda);

#endif