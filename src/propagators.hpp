#ifndef PROPAGATORS_H
#define PROPAGATORS_H

#include <iostream>
#include "armadillo"
#include "operators.hpp"
#include "constants.hpp"

// note: we provide real and complex valued implementations

int get_E_F(arma::sp_mat *E_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_E_F(arma::sp_cx_mat *E_F, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_E_FCF(arma::sp_mat *E_FCF, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_E_FCF(arma::sp_cx_mat *E_FCF, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_F_F(arma::sp_mat *F_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_F_F(arma::sp_cx_mat *F_F, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_F_FCF(arma::sp_mat *F_FCF, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_F_FCF(arma::sp_cx_mat *F_FCF, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel);

int get_R_F(arma::sp_mat *R_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel, int cycle);

int get_R_F(arma::sp_cx_mat *R_F, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel, int cycle);

int get_R_FCF(arma::sp_mat *R_FCF, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel, int cycle);

int get_R_FCF(arma::sp_cx_mat *R_FCF, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel, int cycle);

#endif
