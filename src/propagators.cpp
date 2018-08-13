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
        (*ptrRI[i]) = get_RIl(Nl(i), ml(i));
        (*ptrS[i])  = get_Sl(Nl(i), ml(i));
    }
    ptrA[numLevels-1]       = new arma::sp_mat();
    (*ptrA[numLevels-1])    = get_Al(lambda(numLevels-1), Nl(numLevels-1));
    // compute inverse on coarsest grid
    arma::sp_mat  invA1 = arma::sp_mat(arma::mat((*ptrA[numLevels-1])).i());
    // compute error propagator for F-relaxation
    arma::sp_mat    E_F = (*ptrP[0]) * (*ptrRI[0])
        - (*ptrP[0]) * invA1 * (*ptrRI[0]) * (*ptrA[0]) * (*ptrP[0]) * (*ptrRI[0]);
//    arma::sp_mat    E_F = arma::speye(Nl(1), Nl(1)) - invA1 * (*ptrRI[0]) * (*ptrA[0]) * (*ptrP[0]);
    return E_F;
}

