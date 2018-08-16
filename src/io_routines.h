#ifndef IO_ROUTINES_H
#define IO_ROUTINES_H

#include <iostream>
#include "armadillo"

void export_matrix(arma::sp_mat *m, const std::string& filename, arma::file_type type = arma::raw_ascii);

void export_matrix(arma::sp_cx_mat *m, const std::string& filename, arma::file_type type = arma::raw_ascii);

#endif