#include "error_bound_routines.hpp"

/**
 *  compute bound for error propagator for multilevel MGRIT algorithm
 */
void get_error_propagator_bound(const int bound,                ///< requested bound, see constants.hpp
                                int theoryLevel,                ///< time grid level, where bound is computed
                                const int relax,                ///< relaxation scheme, see constants.hpp
                                Col<int> numberOfTimeSteps,     ///< number of time steps on all time grids
                                Col<int> coarseningFactors,     ///< temporal coarsening factors for levels 0-->1, 1-->2, etc.
                                Col<cx_double> **lambda,        ///< eigenvalues on each time grid
                                Col<double> *&estimate          ///< on return, the estimate for all eigenvalues
                                ){
    int errCode = 0;
    bool pinvSuccess = true;
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
    // setup function pointer for error propagator
    int (*get_E)(arma::sp_cx_mat *, arma::Col<arma::cx_double>, arma::Col<int>, arma::Col<int>, int);
    switch(relax){
        case mgritestimate::F_relaxation:{
            get_E = get_E_F;
            break;
        }
        case mgritestimate::FCF_relaxation:{
            get_E = get_E_FCF;
            break;
        }
        default:{
            cout << ">>>ERROR: Only F- and FCF-relaxation are implemented" << endl;
            throw;
        }
    }
    // get local eigenvalue index range on rank
    if(numberOfSamples < world_size){
        cout << ">>>ERROR: Sample size is too small for number of MPI processes" << endl;
        throw;
    }
    // get index range for local MPI process
    int samplesRankStartIdx;
    int samplesRankStopIdx;
    get_samples_index_range(numberOfSamples, samplesRankStartIdx, samplesRankStopIdx);
    // get estimate
    switch(bound){
        case mgritestimate::upper_bound:{
            errCode = 0;
            for(int evalIdx = samplesRankStartIdx; evalIdx <= samplesRankStopIdx; evalIdx++){
                // for a given spatial mode, get eigenvalues for all levels
                Col<cx_double> lambda_k(numberOfLevels);
                for(int level = 0; level < numberOfLevels; level++){
                    lambda_k(level) = (*lambda[level])(evalIdx);
                }
                sp_cx_mat *E = new sp_cx_mat();
                // get error propagator
                errCode = get_E(E, lambda_k, numberOfTimeSteps, coarseningFactors, 1);
                // compute bound only if time stepper is stable, i.e., \f$\lambda_k < 1\f$
                if(errCode == -1){
                    (*estimate)(evalIdx) = -1;
                    continue;
                }
                (*estimate)(evalIdx) = norm(cx_mat(*E), 2);
            }
            break;
        }
        case mgritestimate::lower_bound:{
            errCode = 0;
            for(int evalIdx = samplesRankStartIdx; evalIdx <= samplesRankStopIdx; evalIdx++){
                // for a given spatial mode, get eigenvalues for all levels
                Col<cx_double> lambda_k(numberOfLevels);
                for(int level = 0; level < numberOfLevels; level++){
                    lambda_k(level) = (*lambda[level])(evalIdx);
                }
                sp_cx_mat *E = new sp_cx_mat();
                cx_mat pinvE;
                pinvSuccess = true;
                // get error propagator
                errCode = get_E(E, lambda_k, numberOfTimeSteps, coarseningFactors, 1);
                // compute bound only if time stepper is stable, i.e., \f$\lambda_k < 1\f$
                if(errCode == -1){
                    (*estimate)(evalIdx) = -1;
                    continue;
                }
                // compute pseudo-inverse of error propagator
                pinvSuccess = pinv(pinvE, cx_mat(*E));
                if(!pinvSuccess){
                    cout << ">>>ERROR: Computing pseudo-inverse of error propagator failed." << endl;
                    throw;
                }
                // l2-norm of pseudo-inverse bounds error propagator from below
                /// \todo should we check for division by zero?
                (*estimate)(evalIdx) = 1.0 / norm(pinvE, 2);
            }
            break;
        }
        case mgritestimate::sqrt_upper_bound:{
            errCode = 0;
            double norm1 = 0.0;
            double normInf = 0.0;
            for(int evalIdx = samplesRankStartIdx; evalIdx <= samplesRankStopIdx; evalIdx++){
                // for a given spatial mode, get eigenvalues for all levels
                Col<cx_double> lambda_k(numberOfLevels);
                for(int level = 0; level < numberOfLevels; level++){
                    lambda_k(level) = (*lambda[level])(evalIdx);
                }
                sp_cx_mat *E = new sp_cx_mat();
                // get error propagator
                errCode = get_E(E, lambda_k, numberOfTimeSteps, coarseningFactors, 1);
                // compute bound only if time stepper is stable, i.e., \f$\lambda_k < 1\f$
                if(errCode == -1){
                    (*estimate)(evalIdx) = -1;
                    continue;
                }
                norm1   = norm(*E, 1);
                normInf = norm(*E, "inf");
                (*estimate)(evalIdx) = sqrt(norm1 * normInf);
            }
            break;
        }
        default:{
            cout << ">>>ERROR: Bound " << bound << " not implemented" << endl;
            throw;
        }
    }
    // communicate data
    communicateBounds(estimate, numberOfSamples, samplesRankStartIdx, samplesRankStopIdx);
}

/**
 *  compute global index range for local MPI rank
 */
void get_samples_index_range(int numberOfSamples,       ///< global number of samples
                             int &samplesRankStartIdx,  ///< start index of samples on MPI rank
                             int &samplesRankStopIdx    ///< stop index of samples on MPI rank
                             ){
    // check MPI environment
    int world_rank;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // compute sample sizes on MPI rank
    int samplesPerRank      = numberOfSamples / world_size;
    int leftover            = numberOfSamples - samplesPerRank * world_size;
    samplesRankStartIdx     = samplesPerRank * world_rank           + leftover * (world_rank > 0);
    samplesRankStopIdx      = samplesPerRank * (world_rank + 1) - 1 + leftover;
    // cout << "Rank " << world_rank << " / " << world_size << ": eigenvalue index range is "
    //      << samplesRankStartIdx << " - " << samplesRankStopIdx << " (" << numberOfSamples << ")"
    //      << endl;
}

/**
 *  collect computed bounds on rank 0
 */
void communicateBounds(Col<double> *&estimate,          ///< on return, the estimate for all eigenvalues
                       int numberOfSamples,             ///< global number of samples
                       int samplesRankStartIdx,         ///< start index of samples on MPI rank
                       int samplesRankStopIdx           ///< stop index of samples on MPI rank
                       ){
    // check MPI environment
    int world_rank;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // send data to rank 0
    int idx = 0;
    if(world_rank == 0){
        int samplesPerRank  = numberOfSamples / world_size; // minimum number of samples on each MPI rank
        int receiveSize     = samplesPerRank;               // size of MPI receive
        int receiveStartIdx = samplesRankStartIdx;          // start index of received samples
        int receiveStopIdx  = samplesRankStopIdx;           // stop index of received samples
        for(int rank = 1; rank < world_size; rank++){
            receiveStartIdx = receiveStopIdx + 1;               // increment start index for MPI rank
            receiveStopIdx  = receiveStopIdx + samplesPerRank;  // increment stop index for MPI rank
            double data[samplesPerRank];
            // post blocking receive
            MPI_Recv(data, receiveSize, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
            idx = 0;
            // unpack data
            for(int evalIdx = receiveStartIdx; evalIdx <= receiveStopIdx; evalIdx++){
                (*estimate)(evalIdx) = data[idx++];
            }
        }
    }else{
        int samplesOnRank = samplesRankStopIdx - samplesRankStartIdx + 1;   // size of MPI send
        double data[samplesOnRank];
        idx = 0;
        // pack data
        for(int evalIdx = samplesRankStartIdx; evalIdx <= samplesRankStopIdx; evalIdx++){
            data[idx++] = (*estimate)(evalIdx);
        }
        // post blocking send
        MPI_Send(data, samplesOnRank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}
