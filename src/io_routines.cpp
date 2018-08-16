#include "io_routines.h"

void export_matrix(arma::sp_mat *m, const std::string& filename, arma::file_type type){
    arma::mat(*m).save(filename, type);
}

void export_matrix(arma::sp_cx_mat *m, const std::string& filename, arma::file_type type){
    arma::mat(arma::real(*m)).save(filename+"_real", type);
    arma::mat(arma::imag(*m)).save(filename+"_imag", type);
}