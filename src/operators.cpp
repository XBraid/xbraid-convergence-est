#include "operators.h"

arma::sp_mat get_Al(double lambda, int Nl){
    // set size
    arma::sp_mat Al(Nl, Nl);
    // set ones on diagonal
    for(int i = 0; i < Nl; i++){
        Al(i,i) = 1.0;
    }
    // set -lambda on off-diagonal (one below)
    for(int i = 0; i < Nl-1; i++){
        Al(i+1,i) = -lambda;
    }
    return Al;
}

arma::sp_cx_mat get_Al(arma::cx_double lambda, int Nl){
    // set size
    arma::sp_cx_mat Al(Nl, Nl);
    // set ones on diagonal
    for(int i = 0; i < Nl; i++){
        Al(i,i) = {1.0,0.0};
    }
    // set -lambda on off-diagonal (one below)
    for(int i = 0; i < Nl-1; i++){
        Al(i+1,i) = -lambda;
    }
    return Al;
}

arma::sp_mat get_Rl(double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_mat Rl(Nl1, Nl);
    // get non-zero values
    arma::vec vals = arma::ones<arma::vec>(ml);
    for(int i = 0; i < ml; i++){
        vals(i) = std::pow(lambda, ml-1-i);
    }
    // set Rl
    int j = 0;
    Rl(j,j) = 1.0;
    for(int i = 1; i < Nl-1; i += ml){
        j++;
        Rl(j,arma::span(i,i+ml-1)) = vals.st();
    }
    return Rl;
}

arma::sp_cx_mat get_Rl(arma::cx_double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_cx_mat Rl(Nl1, Nl);
    // get non-zero values
    arma::cx_vec vals = arma::zeros<arma::cx_vec>(ml);
    for(int i = 0; i < ml; i++){
        vals(i) = std::pow(lambda, ml-1-i);
    }
    // set Rl
    int j = 0;
    Rl(j,j) = {1.0,0.0};
    for(int i = 1; i < Nl-1; i += ml){
        j++;
        Rl(j,arma::span(i,i+ml-1)) = vals.st();
    }
    return Rl;
}

arma::sp_mat get_RIl(double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_mat RIl(Nl1, Nl);
    // set Pl
    int j = 0;
    for(int i = 0; i < Nl1; i++){
        RIl(i,j) = 1.0;
        j += ml;
    }
    return RIl;
}

arma::sp_cx_mat get_RIl(arma::cx_double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_cx_mat RIl(Nl1, Nl);
    // set Pl
    int j = 0;
    for(int i = 0; i < Nl1; i++){
        RIl(i,j) = {1.0,0.0};
        j += ml;
    }
    return RIl;
}

arma::sp_mat get_Sl(double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_mat Sl(Nl, Nl-Nl1);
    // set Sl
    int idx = 0;
    for(int i = 0; i < Nl; i++){
        if((i % ml) > 0){
            Sl(i,idx) = 1.0;
            idx++;
        }
    }
    return Sl;
}

arma::sp_cx_mat get_Sl(arma::cx_double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_cx_mat Sl(Nl, Nl-Nl1);
    // set Sl
    int idx = 0;
    for(int i = 0; i < Nl; i++){
        if((i % ml) > 0){
            Sl(i,idx) = {1.0,0.0};
            idx++;
        }
    }
    return Sl;
}

arma::sp_mat get_Pl(double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_mat Pl(Nl, Nl1);
    // get non-zero values
    arma::vec vals = arma::ones<arma::vec>(ml);
    for(int i = 0; i < ml; i++){
        vals(i) = std::pow(lambda, i);
    }
    // set Pl
    int j = 0;
    for(int i = 0; i < Nl-1; i += ml){
        Pl(arma::span(i,i+ml-1),j) = vals;
        j++;
    }
    Pl(Nl-1,Nl1-1) = 1.0;
    return Pl;
}

arma::sp_cx_mat get_Pl(arma::cx_double lambda, int Nl, int ml){
    // compute coarse-grid size
    int Nl1 = (Nl - 1) / ml + 1;
    // set size
    arma::sp_cx_mat Pl(Nl, Nl1);
    // get non-zero values
    arma::cx_vec vals = arma::ones<arma::cx_vec>(ml);
    for(int i = 0; i < ml; i++){
        vals(i) = std::pow(lambda, i);
    }
    // set Pl
    int j = 0;
    for(int i = 0; i < Nl-1; i += ml){
        Pl(arma::span(i,i+ml-1),j) = vals;
        j++;
    }
    Pl(Nl-1,Nl1-1) = {1.0,0.0};
    return Pl;
}
