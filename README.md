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

To install bootstrap and armadillo, we recommend using SPACK, since SPACK will install all dependencies (blas, lapack, etc.).

## Makefile targets

    default : builds the default binary
    doxygen : generates the Doxygen source code documentation
    load-packages : loads the armadillo library using SPACK
    unload-packages : unloads the armadillo library using SPACK

## Source code documentation

Doxygen is used to automatically generate source code documentation.
See doc/doxygen/html/index.html