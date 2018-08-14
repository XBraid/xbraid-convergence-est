#include "propagators.h"

void get_E_F(arma::sp_mat **E_F, arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml){
    E_F[0] = new arma::sp_mat();
    // get number of levels
    int                 numLevels   = lambda.n_elem;
    // pointers to all operators
    arma::sp_mat *ptrA[numLevels];
    arma::sp_mat *ptrP[numLevels-1];
    arma::sp_mat *ptrR[numLevels-1];
    arma::sp_mat *ptrRI[numLevels-1];
    arma::sp_mat *ptrS[numLevels-1];
    // compute operators
    get_operators(ptrA, ptrR, ptrRI, ptrS, ptrP, lambda, Nl, ml);
    // compute error propagator for F-relaxation, using form: I - term * R0A0P0 - sum * R0A0P0
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
    for(int level = 0; level < numLevels-2; level++){
        term2 = (*ptrS[level+1]).st() * (*ptrA[level+1]) * (*ptrS[level+1]);
        for(int prepostIdx = level; prepostIdx > 1; prepostIdx--){
            term2 = (*ptrP[prepostIdx]) * arma::sp_mat(arma::mat(term2).i()) * (*ptrR[prepostIdx]);
        }
        summ += term2;
    }
    // compute error propagator for F-relaxation
    (*E_F[0]) = arma::speye(Nl(1), Nl(1)) - term * R0A0P0 - summ * R0A0P0;
//    (*E_F) = (*ptrP[0]) * (arma::speye(Nl(1), Nl(1)) - term * R0A0P0 - summ * R0A0P0) * (*ptrRI[0]);
}

