#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"

struct appStruct{
    int world_size;                                                         ///< MPI world size
    int world_rank;                                                         ///< MPI rank ID
    int numberOfTimeSteps_l0;                                               ///< number of time steps on fine grid
    int numberOfLevels;                                                     ///< number of time grid levels
    bool sampleComplexPlane;                                                ///< sample complex plane for \f$z = \delta_t \cdot \eta\f$?
    arma::Col<int> coarseningFactors;                                       ///< temporal coarsening factors between levels 0-->1, 1-->2, etc.
    arma::Col<int> numberOfTimeSteps;                                       ///< number of time steps on levels 0, 1, etc.
    int method;                                                             ///< Runge-Kutta method if spatial eigenvalues are supplied, see constants.hpp
    double min_dteta_real_l0;                                               ///< minimum x-range if sampling spatial eigenvalues
    double max_dteta_real_l0;                                               ///< maximum x-range if sampling spatial eigenvalues
    double min_dteta_imag_l0;                                               ///< minimum y-range if sampling spatial eigenvalues
    double max_dteta_imag_l0;                                               ///< maximum y-range if sampling spatial eigenvalues
    int numberOfRealSamples;                                                ///< number of samples along real axis if sampling spatial eigenvalues
    int numberOfImagSamples;                                                ///< number of samples along complex axis if sampling spatial eigenvalues
    int bound;                                                              ///< bound to evaluate, see constants.hpp
    int theoryLevel;                                                        ///< bound is evaluated on this level
    int relax;                                                              ///< relaxation scheme, see constants.hpp
    bool fileSpatialEigenvalues;                                            ///< indicates whether we are reading spatial eigenvalues from file
    bool filePhiEigenvalues;                                                ///< indicates whether we are reading eigenvalues of \f$Phi_l\f$ from file
    bool fileComplexEigenvalues;                                            ///< indicates whether imported eigenvalues are complex
    std::string fileNameSpatialEigenvaluesReal;                             ///< name of file that contains spatial eigenvalues (real part)
    std::string fileNameSpatialEigenvaluesImag;                             ///< name of file that contains spatial eigenvalues (complex part)
    std::string fileNamePhiEigenvaluesReal;                                 ///< name of file that contains eigenvalues of \f$Phi_l\f$ (real part)
    std::string fileNamePhiEigenvaluesImag;                                 ///< name of file that contains eigenvalues of \f$Phi_l\f$ (complex part)
    arma::Col<double> **dtetar;                                             ///< pointer to spatial eigenvalues at different grid levels (real case)
    arma::Col<double> **lambdar;                                            ///< pointer to eigenvalues of \f$\Phi_l\f$ at different grid levels (complex case)
    arma::Col<arma::cx_double> **dtetac;                                    ///< pointer to spatial eigenvalues at different grid levels (real case)
    arma::Col<arma::cx_double> **lambdac;                                   ///< pointer to eigenvalues of \f$\Phi_l\f$ at different grid levels (complex case)
    arma::Col<double> *estimate;                                            ///< computed estimates
    bool userOutputFile;                                                    ///< indicates whether user has defined output file name
    std::string userOutputFileName;                                         ///< user-defined output file name
};

void initalize_app_struct(appStruct &app);

#endif