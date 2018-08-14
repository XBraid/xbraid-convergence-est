#include "propagators.h"

arma::sp_mat get_E_F(arma::Col<double> lambda, arma::Col<int> Nl, arma::Col<int> ml){
    // get number of levels
    int                 numLevels   = lambda.n_elem;
    std::cout << "Number of levels: " << numLevels << std::endl;
    // pointers to all operators
    arma::sp_mat *ptrA[numLevels];
    arma::sp_mat *ptrP[numLevels-1];
    arma::sp_mat *ptrR[numLevels-1];
    arma::sp_mat *ptrRI[numLevels-1];
    arma::sp_mat *ptrS[numLevels-1];
    // compute operators
    for(int i = 0; i < numLevels-1; i++){
        ptrA[i]     = new arma::sp_mat();
        ptrP[i]     = new arma::sp_mat();
        ptrR[i]     = new arma::sp_mat();
        ptrRI[i]    = new arma::sp_mat();
        ptrS[i]     = new arma::sp_mat();
        (*ptrA[i])  = get_Al(lambda(i), Nl(i));
        (*ptrP[i])  = get_Pl(lambda(i), Nl(i), ml(i));
        (*ptrR[i])  = get_Rl(lambda(i), Nl(i), ml(i));
        (*ptrRI[i]) = get_RIl(lambda(i), Nl(i), ml(i));
        (*ptrS[i])  = get_Sl(lambda(i), Nl(i), ml(i));
    }
    ptrA[numLevels-1]       = new arma::sp_mat();
    (*ptrA[numLevels-1])    = get_Al(lambda(numLevels-1), Nl(numLevels-1));
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
    arma::sp_mat E_F = arma::speye(Nl(1), Nl(1)) - term * R0A0P0 - summ * R0A0P0;
//    E_F = (*ptrP[0]) * E_F * (*ptrRI[0]);
    return E_F;
}

