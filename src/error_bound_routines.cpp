#include "error_bound_routines.hpp"

void get_error_propagator_bound(const int bound,                ///< requested bound, see constants.hpp
                                const int relax,                ///< relaxation scheme, see constants.hpp
                                Col<int> numberOfTimeSteps,     ///< number of time steps on all time grids
                                Col<int> coarseningFactors,     ///< temporal coarsening factors for levels 0-->1, 1-->2, etc.
                                Col<cx_double> **lambda,        ///< eigenvalues on each time grid
                                Col<double> *&estimate          ///< on return, the estimate for all eigenvalues
                                ){
    int errCode = 0;
    // check MPI environment
    int world_rank;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // get number of eigenvalues and levels
    int numberOfSamples = (*lambda[0]).n_elem;
    int numberOfLevels  = numberOfTimeSteps.n_elem;
    estimate = new Col<double>;
    (*estimate).set_size(numberOfSamples);
    // get local eigenvalue index range on rank
    if(numberOfSamples < world_size){
        cout << ">>>ERROR: Sample size is too small for number of MPI processes" << endl;
        throw;
    }
    int samplesPerRank      = numberOfSamples / world_size;
    int leftover            = numberOfSamples - samplesPerRank * world_size;
    int samplesRankStartIdx = samplesPerRank * world_rank           + leftover * (world_rank > 0);
    int samplesRankStopIdx  = samplesPerRank * (world_rank + 1) - 1 + leftover;
    // cout << "Rank " << world_rank << " / " << world_size << ": eigenvalue index range is "
    //      << samplesRankStartIdx << " - " << samplesRankStopIdx << " (" << numberOfSamples << ")"
    //      << endl;
    int samplesOnRank = samplesPerRank + leftover * (world_rank == 0);
    // get estimate
    switch(bound){
        case mgritestimate::upper_bound:{
            switch(relax){
                case mgritestimate::F_relaxation:{
                    errCode = 0;
                    for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
                        // for a given spatial mode, get eigenvalues for all levels
                        Col<cx_double> lambda_k(numberOfLevels);
                        for(int level = 0; level < numberOfLevels; level++){
                            lambda_k(level) = (*lambda[level])(evalIdx);
                        }
                        sp_cx_mat *E_F = new sp_cx_mat();
                        errCode = get_E_F(E_F, lambda_k, numberOfTimeSteps, coarseningFactors);
                        if(errCode == -1){
                            (*estimate)(evalIdx) = -1;
                        }
                        else{
                            (*estimate)(evalIdx) = norm(cx_mat(*E_F), 2);
                        }
                    }
                    break;
                }
                case mgritestimate::FCF_relaxation:{
                    errCode = 0;
                    for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
                        // for a given spatial mode, get eigenvalues for all levels
                        Col<cx_double> lambda_k(numberOfLevels);
                        for(int level = 0; level < numberOfLevels; level++){
                            lambda_k(level) = (*lambda[level])(evalIdx);
                        }
                        sp_cx_mat *E_FCF = new sp_cx_mat();
                        errCode = get_E_FCF(E_FCF, lambda_k, numberOfTimeSteps, coarseningFactors);
                        if(errCode == -1){
                            (*estimate)(evalIdx) = -1;
                        }
                        else{
                            (*estimate)(evalIdx) = norm(cx_mat(*E_FCF), 2);
                        }
                    }
                    break;
                }
                default:{
                    cout << ">>>ERROR: Only F- and FCF-relaxation are implemented" << endl;
                    throw;
                }
            }
            break;
        }
        case mgritestimate::lower_bound:{
            cout << ">>>ERROR: Bound " << bound << " not implemented" << endl;
            throw;
        }
        default:{
            cout << ">>>ERROR: Bound " << bound << " not implemented" << endl;
            throw;
        }
    }
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
                (*estimate)(evalIdx) = data[idx++];
            }
        }
    }else{
        double data[samplesOnRank];
        idx = 0;
        for(int evalIdx = samplesRankStartIdx; evalIdx < samplesRankStopIdx + 1; evalIdx++){
            data[idx++] = (*estimate)(evalIdx);
        }
        MPI_Send(data, samplesOnRank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}
