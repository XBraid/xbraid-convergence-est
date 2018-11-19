# mgrit-convergence

## Operating systems

This code was developed for use with Linux operating systems.

## Installation requirements

The code depends on a number of third-party libraries:

* C++ compiler
* MPI (recommended: MPICH2)
* bootstrap
* armadillo
* optional: Doxygen
* optional: matplotlib (recommended: version >= 1.5.2)

To install bootstrap and armadillo, we recommend using SPACK, since SPACK will install all dependencies (blas, lapack, etc.).

## Makefile targets

    default : builds the default binary
    doxygen : generates the Doxygen source code documentation
    load-packages : loads the recommended libraries using SPACK
    unload-packages : unloads the recommended libraries using SPACK

## Commandline options

A list of commandline options is displayed when running:

    ./bin/mgrit-convergence --help

## Visualization

A Python3 script for visualization of estimated bounds for spatial eigenvalues in the complex plane can be found in the utils/ folder.
For required commandline arguments, run:

    python3 heatmap.py

## Source code documentation

Doxygen is used to automatically generate source code documentation.
See doc/doxygen/html/index.html
