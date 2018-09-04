#include <mpi.h>
#include <iostream>
#include <ctime>
#include "armadillo"
#include "operators.hpp"
#include "propagators.hpp"
#include "io_routines.hpp"
#include "error_bound_routines.hpp"
#include "types.hpp"

using namespace arma;
using namespace std;


int main(int argc, char** argv){
    MPI_Init(NULL, NULL);
    appStruct app;
    MPI_Comm_size(MPI_COMM_WORLD, &app.world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &app.world_rank);
    clock_t begin,end;
    int errCode = 0;
    string filename;

    cout << "Armadillo version: " << arma_version::as_string() << endl;
    cout << "MPI rank: " << app.world_rank << " / " << app.world_size << endl;

    // MGRIT settings - set defaults and parse commandline options
    initalize_app_struct(app);
    parse_commandline_options(app, argc, argv);
    // read eigenvalues based on commandline options
    setget_eigenvalues(app);

    Col<double> *estimate;
    begin = clock();
    // evaluate error propagator for complex eigenvalue case
    if(app.sampleComplexPlane || app.fileComplexEigenvalues){
        get_error_propagator_bound(app.bound, app.theoryLevel, app.relax, app.numberOfTimeSteps, app.coarseningFactors, app.lambdac, estimate);
    // evaluate error propagator for real eigenvalue case
    }else{
        get_error_propagator_bound(app.bound, app.theoryLevel, app.relax, app.numberOfTimeSteps, app.coarseningFactors, app.lambdar, estimate);
    }
    end = clock();
    cout << "Bound on rank " << app.world_rank << " / " << app.world_size << " - Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << " seconds" << endl;
    if(app.world_rank == 0){
        get_default_filename(app.bound, app.relax, &filename);
        mat(*estimate).save(filename, raw_ascii);
    }

    MPI_Finalize();
}
