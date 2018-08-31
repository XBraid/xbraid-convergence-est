#include <mpi.h>
#include <iostream>
#include <ctime>
#include "armadillo"
#include "operators.hpp"
#include "propagators.hpp"
#include "io_routines.hpp"
#include "rk_routines.hpp"
#include "sampling_routines.hpp"
#include "error_bound_routines.hpp"

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
    string filename;

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

    // decide which bound to compute
    const int bound = mgritestimate::tight_twogrid_upper_bound;
    int theoryLevel = 1;

    // compute estimate - F-relaxation
    Col<double> *estimateF;
    begin = clock();
    get_error_propagator_bound(bound, theoryLevel, mgritestimate::F_relaxation, numberOfTimeSteps, coarseningFactors, lambda, estimateF);
    end = clock();
    cout << "F-relaxation on rank " << world_rank << " / " << world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    if(world_rank == 0){
        mat(join_rows(real(*dteta[0]), imag(*dteta[0]))).save("dteta_l0.txt", raw_ascii);
        if(bound == mgritestimate::upper_bound){
            filename = "upper_bound_E_F.txt";
        }else if(bound == mgritestimate::sqrt_upper_bound){
            filename = "sqrt_upper_bound_E_F.txt";
        }else if(bound == mgritestimate::sqrt_expression_upper_bound){
            filename = "sqrt_expression_upper_bound_E_F.txt";
        }else if(bound == mgritestimate::lower_bound){
            filename = "lower_bound_E_F.txt";
        }else{
            filename = "bound_E_F.txt";
        }
        mat(*estimateF).save(filename, raw_ascii);
    }

    // compute estimate - FCF-relaxation
    Col<double> *estimateFCF;
    begin = clock();
    get_error_propagator_bound(bound, theoryLevel, mgritestimate::FCF_relaxation, numberOfTimeSteps, coarseningFactors, lambda, estimateFCF);
    end = clock();
    cout << "FCF-relaxation on rank " << world_rank << " / " << world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    if(world_rank == 0){
        mat(join_rows(real(*dteta[0]), imag(*dteta[0]))).save("dteta_l0.txt", raw_ascii);
        string filename;
        if(bound == mgritestimate::upper_bound){
            filename = "upper_bound_E_FCF.txt";
        }else if(bound == mgritestimate::sqrt_upper_bound){
            filename = "sqrt_upper_bound_E_FCF.txt";
        }else if(bound == mgritestimate::sqrt_expression_upper_bound){
            filename = "sqrt_expression_upper_bound_E_FCF.txt";
        }else if(bound == mgritestimate::lower_bound){
            filename = "lower_bound_E_FCF.txt";
        }else{
            filename = "bound_E_FCF.txt";
        }
        mat(*estimateFCF).save(filename, raw_ascii);
    }

    MPI_Finalize();
}
