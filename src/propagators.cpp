#include "propagators.h"

// get error propagator for V-cycle with F-relaxation and >= 2 grid levels (real eigenvalues)
void get_E_F(arma::sp_mat *E_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel){
    // check for valid arguments
    if((theoryLevel != 0) && (theoryLevel != 1)){
        std::cout << ">>>ERROR: Error propagator only implemented for level 0 and level 1." << std::endl;
        throw;
    }
    if((lambda.n_elem != Nl.n_elem) || (lambda.n_elem != ml.n_elem+1)){
        std::cout << ">>>ERROR: Error propagator encountered invalid definition of number of levels." << std::endl;
    }
    // get number of levels
    int                 numLevels   = lambda.n_elem;
    // pointers to all operators
    arma::sp_mat *ptrA[numLevels];
    arma::sp_mat *ptrR[numLevels-1];
    arma::sp_mat *ptrRI[numLevels-1];
    arma::sp_mat *ptrS[numLevels-1];
    arma::sp_mat *ptrP[numLevels-1];
    // compute operators
    get_operators(ptrA, ptrR, ptrRI, ptrS, ptrP, lambda, Nl, ml);
    // compute error propagator for F-relaxation, using form: I - term * R0A0P0 - summ * R0A0P0
    // pre-compute R0*A0*P0 = RI0*A0*P0
    arma::sp_mat R0A0P0 = (*ptrRI[0]) * (*ptrA[0]) * (*ptrP[0]);
    // pre-compute coarse-grid contribution
    arma::sp_mat term = arma::sp_mat(arma::mat((*ptrA[numLevels-1])).i());
    for(int prepostIdx = numLevels-2; prepostIdx > 0; prepostIdx--){
        term = (*ptrP[prepostIdx]) * term * (*ptrR[prepostIdx]);
    }
    // pre-compute contributions from intermediate levels
    arma::sp_mat  summ(Nl(1),Nl(1));
    arma::sp_mat term2 = arma::sp_mat();
    arma::sp_mat identity = arma::speye(Nl(1),Nl(1));
    for(int level = 0; level < numLevels-2; level++){
        term2 = (*ptrS[level+1]).st() * (*ptrA[level+1]) * (*ptrS[level+1]);
        for(int prepostIdx = level; prepostIdx > 1; prepostIdx--){
            term2 = (*ptrP[prepostIdx]) * arma::sp_mat(arma::mat(term2).i()) * (*ptrR[prepostIdx]);
        }
        summ += term2;
    }
    // compute error propagator for F-relaxation on level 0
    if(theoryLevel == 1){
        (*E_F) = identity - term * R0A0P0 - summ * R0A0P0;
    // compute error propagator for F-relaxation on level 1
    } else if(theoryLevel == 0){
        (*E_F) = (*ptrP[0]) * (identity - term * R0A0P0 - summ * R0A0P0) * (*ptrRI[0]);
    }
}

// get error propagator for V-cycle with F-relaxation and >= 2 grid levels (complex eigenvalues)
void get_E_F(arma::sp_cx_mat *E_F, arma::Col<arma::cx_double> lambda, arma::Col<int> Nl, arma::Col<int> ml, int theoryLevel){
    // check for valid arguments
    if((theoryLevel != 0) && (theoryLevel != 1)){
        std::cout << ">>>ERROR: Error propagator only implemented for level 0 and level 1." << std::endl;
        throw;
    }
    if((lambda.n_elem != Nl.n_elem) || (lambda.n_elem != ml.n_elem+1)){
        std::cout << ">>>ERROR: Error propagator encountered invalid definition of number of levels." << std::endl;
    }
    // get number of levels
    int                 numLevels   = lambda.n_elem;
    // pointers to all operators
    arma::sp_cx_mat *ptrA[numLevels];
    arma::sp_cx_mat *ptrR[numLevels-1];
    arma::sp_cx_mat *ptrRI[numLevels-1];
    arma::sp_cx_mat *ptrS[numLevels-1];
    arma::sp_cx_mat *ptrP[numLevels-1];
    // compute operators
    get_operators(ptrA, ptrR, ptrRI, ptrS, ptrP, lambda, Nl, ml);
    // compute error propagator for F-relaxation, using form: I - term * R0A0P0 - summ * R0A0P0
    // pre-compute R0*A0*P0 = RI0*A0*P0
    arma::sp_cx_mat R0A0P0 = (*ptrRI[0]) * (*ptrA[0]) * (*ptrP[0]);
    // pre-compute coarse-grid contribution
    arma::sp_cx_mat term = arma::sp_cx_mat(arma::cx_mat((*ptrA[numLevels-1])).i());
    for(int prepostIdx = numLevels-2; prepostIdx > 0; prepostIdx--){
        term = (*ptrP[prepostIdx]) * term * (*ptrR[prepostIdx]);
    }
    // pre-compute contributions from intermediate levels
    arma::sp_cx_mat  summ(Nl(1),Nl(1));
    arma::sp_cx_mat term2 = arma::sp_cx_mat();
    arma::sp_cx_mat identity = arma::sp_cx_mat(arma::speye(Nl(1),Nl(1)),arma::speye(Nl(1),Nl(1)).zeros());
    for(int level = 0; level < numLevels-2; level++){
        term2 = (*ptrS[level+1]).st() * (*ptrA[level+1]) * (*ptrS[level+1]);
        for(int prepostIdx = level; prepostIdx > 1; prepostIdx--){
            term2 = (*ptrP[prepostIdx]) * arma::sp_cx_mat(arma::cx_mat(term2).i()) * (*ptrR[prepostIdx]);
        }
        summ += term2;
    }
    // compute error propagator for F-relaxation on level 0
    if(theoryLevel == 1){
        (*E_F) = identity - term * R0A0P0 - summ * R0A0P0;
    // compute error propagator for F-relaxation on level 1
    } else if(theoryLevel == 0){
        (*E_F) = (*ptrP[0]) * (identity - term * R0A0P0 - summ * R0A0P0) * (*ptrRI[0]);
    }
}

