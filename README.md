# XBraid-convergence-est

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
    tests : run tests to check consistency of bounds (runtime: approx. 30 minutes)
    tests-small : run small subset of tests to check consistency of bounds (runtime: approx. 3 minutes)

## Source code documentation

Doxygen is used to automatically generate source code documentation.
See doc/doxygen/html/index.html

## Wiki

See the XBraid-convergence-est [Wiki pages](https://github.com/XBraid/xbraid-convergence-est/wiki) for more details.

## License

XBraid-convergence-est is distributed under the terms of both the MIT license and the Apache License (Version 2.0). Users may choose either license, at their option.

All new contributions must be made under both the MIT and Apache-2.0 licenses.

See [LICENSE-MIT](https://github.com/XBraid/xbraid-convergence-est/blob/master/LICENSE-MIT), [LICENSE-APACHE](https://github.com/XBraid/xbraid-convergence-est/blob/master/LICENSE-APACHE) and [COPYRIGHT](https://github.com/XBraid/xbraid-convergence-est/blob/master/COPYRIGHT) for details.
