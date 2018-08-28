#ifndef RK_ROUTINES_H
#define RK_ROUTINES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"

using namespace rkconst;

void stability_function(const int method, arma::Col<double> *z, arma::Col<double> *rz);

void stability_function(const int method, arma::Col<arma::cx_double> *z, arma::Col<arma::cx_double> *rz);

void stability_function(arma::mat A, arma::Col<double> b, arma::Col<double> c, arma::Col<double> *z, arma::Col<double> *rz);

void stability_function(arma::mat A, arma::Col<double> b, arma::Col<double> c, arma::Col<arma::cx_double> *z, arma::Col<arma::cx_double> *rz);

void get_butcher_tableau(const int method, arma::mat &A, arma::Col<double> &b, arma::Col<double> &c);

#endif