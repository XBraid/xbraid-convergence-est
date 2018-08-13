#ifndef OPERATORS_H
#define OPERATORS_H

#include <iostream>
#include "armadillo"

arma::sp_mat get_Al(double lambda, int Nl);

arma::sp_mat get_Rl(double lambda, int Nl, int ml);

arma::sp_mat get_RIl(int Nl, int ml);

arma::sp_mat get_Sl(int Nl, int ml);

arma::sp_mat get_Pl(double lambda, int Nl, int ml);

#endif
