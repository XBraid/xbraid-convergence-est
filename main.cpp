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

    // get local eigenvalue index range on rank
    int samplesPerRank      = numberOfSamples / world_size;
    int leftover            = numberOfSamples - samplesPerRank * world_size;
    int samplesRankStartIdx = samplesPerRank * world_rank           + leftover * (world_rank > 0);
    int samplesRankStopIdx  = samplesPerRank * (world_rank + 1) - 1 + leftover;
    cout << "Rank " << world_rank << " / " << world_size << ": eigenvalue index range is "
         << samplesRankStartIdx << " - " << samplesRankStopIdx << " (" << numberOfSamples << ")"
         << endl;
    int samplesOnRank = samplesPerRank + leftover * (world_rank == 0);

    // compute estimate - F-relaxation
    errCode = 0;
    begin = clock();
    Col<double> estimateF(numberOfSamples);
    for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
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
    cout << "F-relaxation on rank " << world_rank << " / " << world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    // compute estimate - F-relaxation
    errCode = 0;
    begin = clock();
    Col<double> estimateFCF(numberOfSamples);
    for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
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
    cout << "FCF-relaxation on rank " << world_rank << " / " << world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;

    // send data to rank 0
    int idx = 0;
    if(world_rank == 0){
        int receiveSize     = samplesPerRank;
        int receiveStartIdx = samplesRankStartIdx;
        int receiveStopIdx  = samplesRankStopIdx;
        for(int rank = 1; rank < world_size; rank++){
            receiveStartIdx = receiveStopIdx + 1;
            receiveStopIdx  = receiveStopIdx + samplesPerRank;
            double data[samplesPerRank];
            MPI_Recv(data, receiveSize, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            idx = 0;
            for(int evalIdx = receiveStartIdx; evalIdx < receiveStopIdx + 1; evalIdx++){
                estimateF(evalIdx) = data[idx++];
            }
            MPI_Recv(data, receiveSize, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            idx = 0;
            for(int evalIdx = receiveStartIdx; evalIdx < receiveStopIdx + 1; evalIdx++){
                estimateFCF(evalIdx) = data[idx++];
            }
        }
        estimateF.print();
    }else{
        double data[samplesOnRank];
        idx = 0;
        for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
            data[idx++] = estimateF(evalIdx);
        }
        MPI_Send(data, samplesOnRank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        idx = 0;
        for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
            data[idx++] = estimateFCF(evalIdx);
        }
        MPI_Send(data, samplesOnRank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

//     // save to disk
//     export_matrix(E_F, "E_F", raw_ascii);
//     export_matrix(E_FCF, "E_FCF", raw_ascii);

    MPI_Finalize();
}
