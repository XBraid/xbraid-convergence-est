#ifndef IO_ROUTINES_H
#define IO_ROUTINES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"
#include "types.hpp"
#include "sampling_routines.hpp"
#include "rk_routines.hpp"

using namespace std;

void export_matrix(arma::sp_mat *m, const string& filename, arma::file_type type = arma::raw_ascii);

void export_matrix(arma::sp_cx_mat *m, const string& filename, arma::file_type type = arma::raw_ascii);

void get_default_filename(const int bound, const int relax, string *filename);

int parse_commandline_options(appStruct &app, int argc, char** argv);

void setget_eigenvalues(appStruct &app);

#endif