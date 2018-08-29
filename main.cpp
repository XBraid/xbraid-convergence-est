#include <mpi.h>
#include <iostream>
#include <ctime>
#include "armadillo"
#include "operators.hpp"
#include "propagators.hpp"
#include "io_routines.hpp"
#include "rk_routines.hpp"
#include "sampling_routines.hpp"

using namespace arma;
using namespace std;


int main(int argc, char** argv){
    int world_size;
    int world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    clock_t begin,end;
    int errCode = 0;

    cout << "Armadillo version: " << arma_version::as_string() << endl;
    cout << "MPI rank: " << world_rank << " / " << world_size << endl;

    // MGRIT settings
    int numberOfTimeSteps_l0 = 1025;
    int numberOfLevels = 3;
    Col<int> coarseningFactors;
    Col<int> numberOfTimeSteps;
    coarseningFactors.set_size(numberOfLevels-1);
    numberOfTimeSteps.set_size(numberOfLevels);
    coarseningFactors = {2, 4};
    numberOfTimeSteps(0) = numberOfTimeSteps_l0;
    for(int level = 1; level < numberOfLevels; level++){ numberOfTimeSteps(level) = (numberOfTimeSteps(level-1) - 1) / coarseningFactors(level-1) + 1; }
    // sample complex plane for dt*eta
    const int method = rkconst::A_stable_SDIRK2;
    double min_dteta_real_l0 = -10.0;
    double max_dteta_real_l0 =   1.0;
    double min_dteta_imag_l0 =  -4.0;
    double max_dteta_imag_l0 =   4.0;
    int numberOfRealSamples = 11;
    int numberOfImagSamples = 9;
    int numberOfSamples = numberOfRealSamples * numberOfImagSamples;
    Col<cx_double> *dteta[numberOfLevels];
    Col<cx_double> *lambda[numberOfLevels];
    int multiplyBy = 1;
    for(int level = 0; level < numberOfLevels; level++){
        if(level > 0){
            multiplyBy *= coarseningFactors(level-1);
        }
        sample_complex_plane(method,
                             min_dteta_real_l0*multiplyBy, max_dteta_real_l0*multiplyBy, min_dteta_imag_l0*multiplyBy, max_dteta_imag_l0*multiplyBy,
                             numberOfRealSamples, numberOfImagSamples,
                             dteta[level], lambda[level]);
    }
    // compute estimate - F-relaxation
    errCode = 0;
    cout << "Estimate for F-relaxation" << endl;
    begin = clock();
    Col<double> estimateF(numberOfSamples);
    for(int evalIdx = 0; evalIdx < numberOfSamples; evalIdx++){
        // for a given spatial mode, get eigenvalues for all levels
        Col<cx_double> lambda_k(numberOfLevels);
        for(int level = 0; level < numberOfLevels; level++){
            lambda_k(level) = (*lambda[level])(evalIdx);
        }
        sp_cx_mat *E_F = new sp_cx_mat();
        errCode = get_E_F(E_F, lambda_k, numberOfTimeSteps, coarseningFactors);
        if(errCode == -1){
            estimateF(evalIdx) = -1;
        }
        else{
            estimateF(evalIdx) = norm(cx_mat(*E_F), 2);
        }
    }
    end = clock();
    cout << "Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    // compute estimate - F-relaxation
    errCode = 0;
    cout << "Estimate for FCF-relaxation" << endl;
    begin = clock();
    Col<double> estimateFCF(numberOfSamples);
    for(int evalIdx = 0; evalIdx < numberOfSamples; evalIdx++){
        // for a given spatial mode, get eigenvalues for all levels
        Col<cx_double> lambda_k(numberOfLevels);
        for(int level = 0; level < numberOfLevels; level++){
            lambda_k(level) = (*lambda[level])(evalIdx);
        }
        sp_cx_mat *E_FCF = new sp_cx_mat();
        errCode = get_E_FCF(E_FCF, lambda_k, numberOfTimeSteps, coarseningFactors);
        if(errCode == -1){
            estimateFCF(evalIdx) = -1;
        }
        else{
            estimateFCF(evalIdx) = norm(cx_mat(*E_FCF), 2);
        }
    }
    end = clock();
    cout << "Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;

//     // save to disk
//     export_matrix(E_F, "E_F", raw_ascii);
//     export_matrix(E_FCF, "E_FCF", raw_ascii);

    MPI_Finalize();
}
