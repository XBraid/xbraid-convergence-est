#ifndef IO_ROUTINES_H
#define IO_ROUTINES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"

void export_matrix(arma::sp_mat *m, const std::string& filename, arma::file_type type = arma::raw_ascii);

void export_matrix(arma::sp_cx_mat *m, const std::string& filename, arma::file_type type = arma::raw_ascii);

void get_default_filename(const int bound, const int relax, std::string *filename);

#endif