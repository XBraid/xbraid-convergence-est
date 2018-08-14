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

void            get_operators(arma::sp_mat **ptrA, arma::sp_mat **ptrR, arma::sp_mat **ptrRI, arma::sp_mat **ptrS, arma::sp_mat **ptrP,
                              arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml);

#endif
