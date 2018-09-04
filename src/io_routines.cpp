#include "io_routines.hpp"

void export_matrix(arma::sp_mat *m, const string& filename, arma::file_type type){
    arma::mat(*m).save(filename, type);
}

// check newer Armadillo versions and arma::arma_ascii format
void export_matrix(arma::sp_cx_mat *m, const string& filename, arma::file_type type){
    arma::mat(arma::real(*m)).save(filename+"_real", type);
    arma::mat(arma::imag(*m)).save(filename+"_imag", type);
}

/**
 *  set default filename depending on bound/relaxation type
 */
void get_default_filename(const int bound, const int relax, string *filename){
    switch(relax){
        case mgritestimate::F_relaxation:{
            if(bound == mgritestimate::upper_bound){
                *filename = "upper_bound_E_F.txt";
            }else if(bound == mgritestimate::sqrt_upper_bound){
                *filename = "sqrt_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::sqrt_expression_upper_bound){
                *filename = "sqrt_expression_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::tight_twogrid_upper_bound){
                *filename = "tight_twogrid_upper_bound_E_F.txt";
            }else if(bound == mgritestimate::lower_bound){
                *filename = "lower_bound_E_F.txt";
            }else{
                *filename = "bound_E_F.txt";
            }
            break;
        }
        case mgritestimate::FCF_relaxation:{
            if(bound == mgritestimate::upper_bound){
                *filename = "upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::sqrt_upper_bound){
                *filename = "sqrt_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::sqrt_expression_upper_bound){
                *filename = "sqrt_expression_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::tight_twogrid_upper_bound){
                *filename = "tight_twogrid_upper_bound_E_FCF.txt";
            }else if(bound == mgritestimate::lower_bound){
                *filename = "lower_bound_E_FCF.txt";
            }else{
                *filename = "bound_E_FCF.txt";
            }
            break;
        }
        default:{
            cout << ">>>ERROR: Unknown relaxation type when setting default filename." << endl;
            throw;
        }
    }
}

/**
 *  Read commandline options (default arguments):
 *
 *      --help                                                  prints all commandline options
 *      --number-of-timesteps 1025                              number of timesteps on fine grid (must contain time point corresponding to initial condition)
 *      --number-of-time-grids 2                                number of levels in time grid hierarchy
 *      --coarsening-factors 2                                  temporal coarsening factors for level 0-->1, 1-->2, ...
 *                                                              note: must not be used before option --number-of-time-grids
 *      --runge-kutta-method L_stable_SDIRK1                    Runge-Kutta method, if spatial eigenvalues are supplied, see constants.hpp
 *      --sample-complex-plane -10.0 1.0 -4.0 4.0               sampling of complex plane with ranges [-10.0, 1.0] x [-4.0, 4.0]
 *      --complex-plane-sample-size 12 9                        sampling of real and complex axis using 12 and 9 points
 *      --file-spatial-real-eigenvalues dteta_real.txt          name of file that contains real part of spatial eigenvalues
 *      --file-spatial-complex-eigenvalues dteta_imag.txt       name of file that contains complex part of spatial eigenvalues
 *      --file-phi-real-eigenvalues lambda_real.txt             name of file that contains real part of eigenvalues of Phi
 *      --file-phi-complex-eigenvalues lambda_imag.txt          name of file that contains complex part of eigenvalues of Phi
 *      --bound tight_twogrid_upper_bound                       bound to evaluate, see constants.hpp
 *      --bound-on-level 1                                      bound is evaluated on this level (options: 0, 1)
 *      --relaxation-scheme F_relaxation                        relaxation scheme, see constants.hpp
 */
int parse_commandline_options(appStruct &app, int argc, char** argv){
    // need to set number of levels before reading temporal coarsening factors
    bool numberOfLevelsSet = false;
    // parse commandline arguments
    for(int argIdx = 1; argIdx < argc; argIdx++){
        if(string(argv[argIdx]) == "--help"){
            if(app.world_rank == 0){
                cout << "Commandline options:" << endl << endl;
                cout << "    --help                                                  prints all commandline options" << endl;
                cout << "    --number-of-timesteps 1025                              number of timesteps on fine grid (must contain time point corresponding to initial condition)" << endl;
                cout << "    --number-of-time-grids 2                                number of levels in time grid hierarchy" << endl;
                cout << "    --coarsening-factors 2                                  temporal coarsening factors for level 0-->1, 1-->2, ..." << endl;
                cout << "                                                            note: must not be used before option --number-of-time-grids" << endl;
                cout << "    --runge-kutta-method L_stable_SDIRK1                    Runge-Kutta method, if spatial eigenvalues are supplied, see constants.hpp" << endl;
                cout << "    --sample-complex-plane -10.0 1.0 -4.0 4.0               sampling of complex plane with ranges [-10.0, 1.0] x [-4.0, 4.0]" << endl;
                cout << "    --complex-plane-sample-size 12 9                        sampling of real and complex axis using 12 and 9 points" << endl;
                cout << "    --file-spatial-real-eigenvalues dteta_real.txt          name of file that contains real part of spatial eigenvalues" << endl;
                cout << "    --file-spatial-complex-eigenvalues dteta_imag.txt       name of file that contains complex part of spatial eigenvalues" << endl;
                cout << "    --file-phi-real-eigenvalues lambda_real.txt             name of file that contains real part of eigenvalues of Phi" << endl;
                cout << "    --file-phi-complex-eigenvalues lambda_imag.txt          name of file that contains complex part of eigenvalues of Phi" << endl;
                cout << "    --bound tight_twogrid_upper_bound                       bound to evaluate, see constants.hpp" << endl;
                cout << "    --bound-on-level 1                                      bound is evaluated on this level (options: 0, 1)" << endl;
                cout << "    --relaxation-scheme F_relaxation                        relaxation scheme, see constants.hpp" << endl;
                cout << endl;
            }
            return 1;
        }else if(string(argv[argIdx]) == "--number-of-timesteps"){
            app.numberOfTimeSteps_l0 = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--number-of-time-grids"){
            app.numberOfLevels = stoi(argv[++argIdx]);
            app.coarseningFactors.set_size(app.numberOfLevels-1);
            app.coarseningFactors.fill(2);
            app.numberOfTimeSteps.set_size(app.numberOfLevels);
            numberOfLevelsSet = true;
        }else if(string(argv[argIdx]) == "--coarsening-factors"){
            if(!numberOfLevelsSet){
                cout << ">>>ERROR: Please set number of time grid levels before specifying coarsening factors." << endl;
                throw;
            }
            for(int level = 0; level < app.numberOfLevels-1; level++){
                app.coarseningFactors(level) = stoi(argv[++argIdx]);
            }
        }else if(string(argv[argIdx]) == "--runge-kutta-method"){
            argIdx++;
            if(string(argv[argIdx]) == "A_stable_SDIRK2"){
                app.method = rkconst::A_stable_SDIRK2;
            }else if(string(argv[argIdx]) == "A_stable_SDIRK3"){
                app.method = rkconst::A_stable_SDIRK4;
            }else if(string(argv[argIdx]) == "A_stable_SDIRK4"){
                app.method = rkconst::A_stable_SDIRK3;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK1"){
                app.method = rkconst::L_stable_SDIRK1;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK2"){
                app.method = rkconst::L_stable_SDIRK2;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK3"){
                app.method = rkconst::L_stable_SDIRK3;
            }else if(string(argv[argIdx]) == "L_stable_SDIRK4"){
                app.method = rkconst::L_stable_SDIRK4;
            }else{
                cout << ">>>ERROR: Unknown Runge-Kutta method " << string(argv[argIdx]) << "." << endl;
                throw;
            }
        }else if(string(argv[argIdx]) == "--sample-complex-plane"){
            app.sampleComplexPlane = true;
            app.min_dteta_real_l0 = stod(argv[++argIdx]);
            app.max_dteta_real_l0 = stod(argv[++argIdx]);
            app.min_dteta_imag_l0 = stod(argv[++argIdx]);
            app.max_dteta_imag_l0 = stod(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--complex-plane-sample-size"){
            app.sampleComplexPlane = true;
            app.numberOfRealSamples = stoi(argv[++argIdx]);
            app.numberOfImagSamples = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-spatial-real-eigenvalues"){
            app.fileSpatialEigenvalues = true;
            app.fileNameSpatialEigenvaluesReal = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-spatial-complex-eigenvalues"){
            app.fileSpatialEigenvalues = true;
            app.fileComplexEigenvalues = true;
            app.fileNameSpatialEigenvaluesReal = string(argv[++argIdx]);
            app.fileNameSpatialEigenvaluesImag = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-phi-real-eigenvalues"){
            app.filePhiEigenvalues = true;
            app.fileNamePhiEigenvaluesReal = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--file-phi-complex-eigenvalues"){
            app.filePhiEigenvalues = true;
            app.fileComplexEigenvalues = true;
            app.fileNamePhiEigenvaluesReal = string(argv[++argIdx]);
            app.fileNamePhiEigenvaluesImag = string(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--bound"){
            argIdx++;
            if(string(argv[argIdx]) == "upper_bound"){
                app.bound = mgritestimate::upper_bound;
            }else if(string(argv[argIdx]) == "sqrt_upper_bound"){
                app.bound = mgritestimate::sqrt_upper_bound;
            }else if(string(argv[argIdx]) == "sqrt_expression_upper_bound"){
                app.bound = mgritestimate::sqrt_expression_upper_bound;
            }else if(string(argv[argIdx]) == "tight_twogrid_upper_bound"){
                app.bound = mgritestimate::tight_twogrid_upper_bound;
            }else if(string(argv[argIdx]) == "lower_bound"){
                app.bound = mgritestimate::lower_bound;
            }else if(string(argv[argIdx]) == "sqrt_lower_bound"){
                app.bound = mgritestimate::sqrt_lower_bound;
            }else if(string(argv[argIdx]) == "tight_twogrid_lower_bound"){
                app.bound = mgritestimate::tight_twogrid_lower_bound;
            }else{
                cout << "ERROR: Unknown bound " << string(argv[argIdx]) << "." << endl;
            }
        }else if(string(argv[argIdx]) == "--bound-on-level"){
            app.theoryLevel = stoi(argv[++argIdx]);
        }else if(string(argv[argIdx]) == "--relaxation-scheme"){
            argIdx++;
            if(string(argv[argIdx]) == "F_relaxation"){
                app.relax = mgritestimate::F_relaxation;
            }else if(string(argv[argIdx]) == "FCF_relaxation"){
                app.relax = mgritestimate::FCF_relaxation;
            }
        }else{
            cout << ">>>ERROR: Unknown argument " << string(argv[argIdx]) << endl;
            throw;
        }
    }
    // sanity checks for user-defined parameters
    if((app.sampleComplexPlane && app.fileSpatialEigenvalues)
        || (app.sampleComplexPlane && app.filePhiEigenvalues)
        || (app.fileSpatialEigenvalues && app.filePhiEigenvalues)){
        cout << ">>>ERROR: Defined multiple input sources." << endl;
        throw;
    }
    // set number of time steps based on fine grid and coarsening factors
    app.numberOfTimeSteps(0) = app.numberOfTimeSteps_l0;
    for(int level = 1; level < app.numberOfLevels; level++){ app.numberOfTimeSteps(level) = (app.numberOfTimeSteps(level-1) - 1) / app.coarseningFactors(level-1) + 1; }
    return 0;
}

/**
 *  computes or imports real/complex eigenvalues
 */
void setget_eigenvalues(appStruct &app){
    // spatial and Phi eigenvalues
    if(app.sampleComplexPlane || app.fileComplexEigenvalues){
        app.dtetac  = new arma::Col<arma::cx_double>*[app.numberOfLevels];
        app.lambdac = new arma::Col<arma::cx_double>*[app.numberOfLevels];
    }else{
        app.dtetar  = new arma::Col<double>*[app.numberOfLevels];
        app.lambdar = new arma::Col<double>*[app.numberOfLevels];
    }
    // sample complex plane for dt*eta
    if(app.sampleComplexPlane){
        int multiplyBy = 1;
        for(int level = 0; level < app.numberOfLevels; level++){
            if(level > 0){
                multiplyBy *= app.coarseningFactors(level-1);
            }
            sample_complex_plane(app.method,
                                 app.min_dteta_real_l0*multiplyBy, app.max_dteta_real_l0*multiplyBy, app.min_dteta_imag_l0*multiplyBy, app.max_dteta_imag_l0*multiplyBy,
                                 app.numberOfRealSamples, app.numberOfImagSamples,
                                 app.dtetac[level], app.lambdac[level]);
        }
        arma::mat(join_rows(real(*app.dtetac[0]), imag(*app.dtetac[0]))).save("dteta_l0.txt", arma::raw_ascii);
    // import spatial eigenvalues
    }else if(app.fileSpatialEigenvalues){
        if(app.fileComplexEigenvalues){
            // read real and complex part of spatial eigenvalues
            arma::mat realEigs;
            arma::mat imagEigs;
            realEigs.load(app.fileNameSpatialEigenvaluesReal);
            imagEigs.load(app.fileNameSpatialEigenvaluesImag);
            // check whether we have at least as many columns as number of levels
            if(realEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            if(imagEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < app.numberOfLevels; level++){
                arma::Col<arma::cx_double> dtetac(realEigs.col(level), imagEigs.col(level));
                (*app.dtetac)[level] = dtetac;
                stability_function(app.method, app.dtetac[level], app.lambdac[level]);
            }
        }else{
            // read real part of spatial eigenvalues
            arma::mat realEigs;
            realEigs.load(app.fileNameSpatialEigenvaluesReal);
            // check whether we have at least as many columns as number of levels
            if(realEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            // read/store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < app.numberOfLevels; level++){
                arma::Col<double> dteta(realEigs.col(level));
                (*app.dtetar)[level] = dteta;
                stability_function(app.method, app.dtetar[level], app.lambdar[level]);
            }
        }
    // import Phi eigenvalues
    }else if(app.filePhiEigenvalues){
        if(app.fileComplexEigenvalues){
            // read real and complex part of spatial eigenvalues
            arma::mat realEigs;
            arma::mat imagEigs;
            realEigs.load(app.fileNamePhiEigenvaluesReal);
            imagEigs.load(app.fileNamePhiEigenvaluesImag);
            // check whether we have at least as many columns as number of levels
            if(realEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            if(imagEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < app.numberOfLevels; level++){
                arma::Col<arma::cx_double> lambda(realEigs.col(level), imagEigs.col(level));
                (*app.lambdac)[level] = lambda;
            }
        }else{
            // read real part of spatial eigenvalues
            arma::mat realEigs;
            realEigs.load(app.fileNamePhiEigenvaluesReal);
            // check whether we have at least as many columns as number of levels
            if(realEigs.n_cols < app.numberOfLevels){
                cout << ">>>ERROR: Supplied eigenvalues are invalid for " << app.numberOfLevels << " levels." << endl;
                throw;
            }
            // store spatial eigenvalues and compute eigenvalues of Phi for given method
            for(int level = 0; level < app.numberOfLevels; level++){
                arma::Col<double> lambda(realEigs.col(level));
                (*app.lambdar)[level] = lambda;
            }
        }
    }
    // if(app.sampleComplexPlane || app.fileComplexEigenvalues){
    //     (*app.estimate).set_size((*app.lambdac[0]).n_cols);
    // }else{
    //     (*app.estimate).set_size((*app.lambdar[0]).n_cols);
    // }
}
