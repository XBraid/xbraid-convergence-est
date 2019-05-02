#!/bin/bash

np=$1
exe=xbraid-convergence-est
# testlambda=stenosed-valve-flow-2d/lambda_nl-3_dt-0.001_nop.txt
testlambda=stenosed-valve-flow-2d/lambda_nl-3_dt-0.001_nop_head.txt
here=$(pwd)
tol=1.0e-12

cd testing/

numTests=0
numTestsPassed=0
numTestsFailed=0

# note: we are using lambda's that correspond to ml0=2, ml1=2
#       but we use them with other values of ml0 / ml1 as well.
#       thus, we cannot expect any reasonable values
#       and in particular not: convergence estimate smaller than 1
for ml0 in 2 4 8
do
    for ml1 in 2 4 8
    do
        for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
        do
            mpirun -np $np $exe --number-of-timesteps       1025                \
                                --number-of-time-grids      3                   \
                                --coarsening-factors        $ml0 $ml1           \
                                --file-phi-real-eigenvalues $testlambda         \
                                --bound                     $bound              \
                                --bound-on-level            1                   \
                                --F-cycle                                       \
                                --relaxation-scheme         F_relaxation        \
                                --output-file               TEST-$bound         &> /dev/null
        done
        val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
        val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
        absdiff=`python -c "print abs($val0-$val1)<$tol"`
        numTests=$((numTests+1))
        if [ "$absdiff" = "True" ]; then
            numTestsPassed=$((numTestsPassed+1))
        else
            echo ">>>INFO: Test failed for parameters: m_0 = $ml0, m_1 = $ml1. $val0 != $val1"
            numTestsFailed=$((numTestsFailed+1))
        fi
    done
done

echo ""
echo ">>>INFO: Passed tests: $numTestsPassed / $numTests"
echo ""

cd $here
