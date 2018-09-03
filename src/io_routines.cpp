#include "io_routines.hpp"

void export_matrix(arma::sp_mat *m, const std::string& filename, arma::file_type type){
    arma::mat(*m).save(filename, type);
}

// check newer Armadillo versions and arma::arma_ascii format
void export_matrix(arma::sp_cx_mat *m, const std::string& filename, arma::file_type type){
    arma::mat(arma::real(*m)).save(filename+"_real", type);
    arma::mat(arma::imag(*m)).save(filename+"_imag", type);
}

/**
 *  set default filename depending on bound/relaxation type
 */
void get_default_filename(const int bound, const int relax, std::string *filename){
    switch(relax){
        case mgritestimate::F_relaxation:{
            if(bound == mgritestimate::upper_bound){
                *filename = "upper_bound_E_F.txt";
            }else if(bound == mgritestimate::sqrt_upper_bound){
                *filename = "sqrt_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::sqrt_expression_upper_bound){
                *filename = "sqrt_expression_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::tight_twogrid_upper_bound){
                *filename = "tight_twogrid_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::lower_bound){
                *filename = "lower_bound_E_F.txt";
            }else{
                *filename = "bound_E_F.txt";
            }
            break;
        }
        case mgritestimate::FCF_relaxation:{
            if(bound == mgritestimate::upper_bound){
                *filename = "upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::sqrt_upper_bound){
                *filename = "sqrt_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::sqrt_expression_upper_bound){
                *filename = "sqrt_expression_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::tight_twogrid_upper_bound){
                *filename = "tight_twogrid_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::lower_bound){
                *filename = "lower_bound_E_FCF.txt";
            }else{
                *filename = "bound_E_FCF.txt";
            }
            break;
        }
        default:{
            std::cout << ">>>ERROR: Unknown relaxation type when setting default filename." << std::endl;
            throw;
        }
    }
}