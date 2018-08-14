#ifndef OPERATORS_H
#define OPERATORS_H

#include <iostream>
#include "armadillo"

arma::sp_mat    get_Al(double lambda, int Nl);

arma::sp_cx_mat get_Al(arma::cx_double lambda, int Nl);

arma::sp_mat    get_Rl(double lambda, int Nl, int ml);

arma::sp_cx_mat get_Rl(arma::cx_double lambda, int Nl, int ml);

arma::sp_mat    get_RIl(double lambda, int Nl, int ml);

arma::sp_cx_mat get_RIl(arma::cx_double lambda, int Nl, int ml);

arma::sp_mat    get_Sl(double lambda, int Nl, int ml);

arma::sp_cx_mat get_Sl(arma::cx_double lambda, int Nl, int ml);

arma::sp_mat    get_Pl(double lambda, int Nl, int ml);

arma::sp_cx_mat get_Pl(arma::cx_double lambda, int Nl, int ml);

#endif
