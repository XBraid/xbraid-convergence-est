#ifndef PROPAGATORS_H
#define PROPAGATORS_H

#include <iostream>
#include "armadillo"
#include "operators.h"

// note: we provide real and complex valued implementations

void get_E_F(arma::sp_mat **E_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml);

#endif
