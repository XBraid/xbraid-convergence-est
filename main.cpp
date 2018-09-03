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

    // MGRIT settings - set defaults and parse commandline
    int numberOfTimeSteps_l0 = 1025;
    int numberOfLevels = 2;
    bool numberOfLevelsSet = false;
    bool sampleComplexPlane = true;
    Col<int> coarseningFactors;
    Col<int> numberOfTimeSteps;
    coarseningFactors.set_size(numberOfLevels-1);
    coarseningFactors.fill(2);
    numberOfTimeSteps.set_size(numberOfLevels);
    int method = rkconst::A_stable_SDIRK2;
    double min_dteta_real_l0 = -10.0;
    double max_dteta_real_l0 =   1.0;
    double min_dteta_imag_l0 =  -4.0;
    double max_dteta_imag_l0 =   4.0;
    int numberOfRealSamples = 11;
    int numberOfImagSamples = 9;
    int bound = mgritestimate::sqrt_expression_upper_bound;
    int theoryLevel = 1;
    int relax = mgritestimate::F_relaxation;
    bool fileSpatialEigenvalues = false;
    bool filePhiEigenvalues = false;
    bool fileComplexEigenvalues = false;
    string fileNameSpatialEigenvaluesReal;
    string fileNameSpatialEigenvaluesImag;
    string fileNamePhiEigenvaluesReal;
    string fileNamePhiEigenvaluesImag;
    for(int argIdx = 1; argIdx < argc; argIdx++){
        cout << argv[argIdx] << endl;
        if(string(argv[argIdx]) == "--number-of-timesteps"){
            numberOfTimeSteps_l0 = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--number-of-time-grids"){
            numberOfLevels = stoi(argv[++argIdx]);
            coarseningFactors.set_size(numberOfLevels-1);
            coarseningFactors.fill(2);
            numberOfTimeSteps.set_size(numberOfLevels);
            numberOfLevelsSet = true;
        }else if(string(argv[argIdx]) == "--coarsening-factors"){
            if(!numberOfLevelsSet){
                cout << ">>>ERROR: Please set number of time grid levels before specifying coarsening factors." << endl;
                throw;
            }
            for(int level = 0; level < numberOfLevels-1; level++){
                coarseningFactors(level) = stoi(argv[++argIdx]);
            }
        }else if(string(argv[argIdx]) == "--runge-kutta-method"){
            argIdx++;
            if(string(argv[argIdx]) == "A_stable_SDIRK2"){
                method = rkconst::A_stable_SDIRK2;
            }else if(string(argv[argIdx]) == "A_stable_SDIRK3"){
                method = rkconst::A_stable_SDIRK4;
            }else if(string(argv[argIdx]) == "A_stable_SDIRK4"){
                method = rkconst::A_stable_SDIRK3;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK1"){
                method = rkconst::L_stable_SDIRK1;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK2"){
                method = rkconst::L_stable_SDIRK2;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK3"){
                method = rkconst::L_stable_SDIRK3;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK4"){
                method = rkconst::L_stable_SDIRK4;
            }else{
                cout << ">>>ERROR: Unknown Runge-Kutta method " << string(argv[argIdx]) << "." << endl;
                throw;
            }
        }else if(string(argv[argIdx]) == "--sample-complex-plane"){
            sampleComplexPlane = true;
            min_dteta_real_l0 = stod(argv[++argIdx]);
            max_dteta_real_l0 = stod(argv[++argIdx]);
            min_dteta_imag_l0 = stod(argv[++argIdx]);
            max_dteta_imag_l0 = stod(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--complex-plane-sample-number"){
            sampleComplexPlane = true;
            numberOfRealSamples = stoi(argv[++argIdx]);
            numberOfImagSamples = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-spatial-real-eigenvalues"){
            fileSpatialEigenvalues = true;
            fileNameSpatialEigenvaluesReal = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-spatial-complex-eigenvalues"){
            fileSpatialEigenvalues = true;
            fileComplexEigenvalues = true;
            fileNameSpatialEigenvaluesReal = string(argv[++argIdx]);
            fileNameSpatialEigenvaluesImag = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-phi-real-eigenvalues"){
            filePhiEigenvalues = true;
            fileNamePhiEigenvaluesReal = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-phi-complex-eigenvalues"){
            filePhiEigenvalues = true;
            fileComplexEigenvalues = true;
            fileNamePhiEigenvaluesReal = string(argv[++argIdx]);
            fileNamePhiEigenvaluesImag = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--bound"){
            argIdx++;
            if(string(argv[argIdx]) == "upper_bound"){
                bound = mgritestimate::upper_bound;
            }else if(string(argv[argIdx]) == "sqrt_upper_bound"){
                bound = mgritestimate::sqrt_upper_bound;
            }else if(string(argv[argIdx]) == "sqrt_expression_upper_bound"){
                bound = mgritestimate::sqrt_expression_upper_bound;
            }else if(string(argv[argIdx]) == "tight_twogrid_upper_bound"){
                bound = mgritestimate::tight_twogrid_upper_bound;
            }else if(string(argv[argIdx]) == "lower_bound"){
                bound = mgritestimate::lower_bound;
            }else if(string(argv[argIdx]) == "sqrt_lower_bound"){
                bound = mgritestimate::sqrt_lower_bound;
            }else if(string(argv[argIdx]) == "tight_twogrid_lower_bound"){
                bound = mgritestimate::tight_twogrid_lower_bound;
            }else{
                cout << "ERROR: Unknown bound " << string(argv[argIdx]) << "." << endl;
            }
        }else if(string(argv[argIdx]) == "--bound-on-level"){
            theoryLevel = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--relaxation-scheme"){
            argIdx++;
            if(string(argv[argIdx]) == "F_relaxation"){
                relax = mgritestimate::F_relaxation;
            }else if(string(argv[argIdx]) == "FCF_relaxation"){
                relax = mgritestimate::FCF_relaxation;
            }
        }else{
            cout << ">>>ERROR: Unknown argument " << string(argv[argIdx]) << endl;
            throw;
        }
    }
    // sanity checks for user-defined parameters
    if((sampleComplexPlane && fileSpatialEigenvalues)
        || (sampleComplexPlane && filePhiEigenvalues)
        || (fileSpatialEigenvalues && filePhiEigenvalues)){
        cout << ">>>ERROR: Defined multiple input sources." << endl;
        throw;
    }
    // set number of time steps based on fine grid and coarsening factors
    numberOfTimeSteps(0) = numberOfTimeSteps_l0;
    for(int level = 1; level < numberOfLevels; level++){ numberOfTimeSteps(level) = (numberOfTimeSteps(level-1) - 1) / coarseningFactors(level-1) + 1; }
    // sample complex plane for dt*eta
    /// \todo conditional typedef for efficiency
    Col<cx_double> *dteta[numberOfLevels];
    Col<cx_double> *lambda[numberOfLevels];
    if(sampleComplexPlane){
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
        mat(join_rows(real(*dteta[0]), imag(*dteta[0]))).save("dteta_l0.txt", raw_ascii);
    }else if(fileSpatialEigenvalues){
        if(fileComplexEigenvalues){
            // read real and complex part of spatial eigenvalues
            mat realEigs;
            mat imagEigs;
            realEigs.load(fileNameSpatialEigenvaluesReal);
            imagEigs.load(fileNameSpatialEigenvaluesImag);
            // check whether we have at least as many columns as number of levels
            if(realEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            if(imagEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < numberOfLevels; level++){
                (*dteta)[level] = Col<cx_double>(realEigs[:, level], imagEigs[:, level]);
                stability_function(method, dteta[level], lambda[level]);
            }
        }else{
            // read real part of spatial eigenvalues
            mat realEigs;
            realEigs.load(fileNameSpatialEigenvaluesReal);
            // check whether we have at least as many columns as number of levels
            if(realEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            // read/store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < numberOfLevels; level++){
                (*dteta)[level] = Col<double>(realEigs[:, level]);
                stability_function(method, dteta[level], lambda[level]);
            }
        }
    }else if(filePhiEigenvalues){
        if(fileComplexEigenvalues){
            // read real and complex part of spatial eigenvalues
            mat realEigs;
            mat imagEigs;
            realEigs.load(fileNamePhiEigenvaluesReal);
            imagEigs.load(fileNamePhiEigenvaluesImag);
            // check whether we have at least as many columns as number of levels
            if(realEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            if(imagEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < numberOfLevels; level++){
                (*lambda)[level] = Col<cx_double>(realEigs[:, level], imagEigs[:, level]);
            }
        }else{
            // read real part of spatial eigenvalues
            mat realEigs;
            realEigs.load(fileNamePhiEigenvaluesReal);
            // check whether we have at least as many columns as number of levels
            if(realEigs.num_cols() < numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << string(numberOfLevels) << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < numberOfLevels; level++){
                (*lambda)[level] = Col<double>(realEigs[:, level]);
            }
        }
    }

    Col<double> *estimate;
    begin = clock();
    get_error_propagator_bound(bound, theoryLevel, relax, numberOfTimeSteps, coarseningFactors, lambda, estimate);
    end = clock();
    cout << "Bound on rank " << world_rank << " / " << world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    if(world_rank == 0){
        get_default_filename(bound, relax, &filename);
        mat(*estimate).save(filename, raw_ascii);
    }

    MPI_Finalize();
}
