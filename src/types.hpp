#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include "armadillo"
#include "constants.hpp"

struct appStruct{
    int world_size;
    int world_rank;
    int numberOfTimeSteps_l0;
    int numberOfLevels;
    bool sampleComplexPlane;
    arma::Col<int> coarseningFactors;
    arma::Col<int> numberOfTimeSteps;
    int method;
    double min_dteta_real_l0;
    double max_dteta_real_l0;
    double min_dteta_imag_l0;
    double max_dteta_imag_l0;
    int numberOfRealSamples;
    int numberOfImagSamples;
    int bound;
    int theoryLevel;
    int relax;
    bool fileSpatialEigenvalues;
    bool filePhiEigenvalues;
    bool fileComplexEigenvalues;
    std::string fileNameSpatialEigenvaluesReal;
    std::string fileNameSpatialEigenvaluesImag;
    std::string fileNamePhiEigenvaluesReal;
    std::string fileNamePhiEigenvaluesImag;
    arma::Col<double> **dtetar;
    arma::Col<double> **lambdar;
    arma::Col<arma::cx_double> **dtetac;
    arma::Col<arma::cx_double> **lambdac;
};

void initalize_app_struct(appStruct &app);

#endif