#!/bin/bash
#
# usage: bash testing/run_tests.sh np {V-cycle} {F-cycle} {full}
#
#           np      - number of processors provided as first argument
#           V-cycle - argument to select V-cycle tests
#           F-cycle - argument to select F-cycle tests
#           full    - optional argument to run larger number of tests
#

start=`date +%s`

full=false
fcyc=false
vcyc=false
ps=( 1 )
if [[ "$*" == *full* ]]
then
    full=true
    ps=( 1 2 3 )
else
    echo ">>>INFO: Running subset of tests."
fi
if [[ "$*" == *F-cycle* ]]
then
    fcyc=true
    echo ">>>INFO: F-cycle tests selected."
fi
if [[ "$*" == *V-cycle* ]]
then
    vcyc=true
    echo ">>>INFO: V-cycle tests selected."
fi
echo ""

np=$1
exe=xbraid-convergence-est
testdir=phiEvals
here=$(pwd)
tol=1.0e-4

cd testing/

numTests=0
numTestsPassed=0
numTestsFailed=0

if [[ $fcyc == true ]]; then
    echo ">>>INFO: Running 3-level F-cycle tests."
    ########################################################################################################################
    # compare 3-level F-cycle with F-relaxation: inequality bound with inequality bound expression
    nl=3
    ########################################################################################################################
    # note: we are using lambda's that correspond to ml0=2, ml1=2
    #       but we use them with other values of ml0 / ml1 as well.
    #       thus, we cannot expect any reasonable values
    #       and in particular not: convergence estimate smaller than 1
    for (( i = 0; i < ${#ps[@]}; ++i ))
    do
        p=${ps[i]}
        # first, test even coarsening factors
        testlambdar=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_real.txt
        testlambdai=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_imag.txt
        for lambdari in "--file-phi-real-eigenvalues $testlambdar" "--file-phi-complex-eigenvalues $testlambdar $testlambdai"
        do
            for ml0 in 2 4 8
            do
                for ml1 in 2 4 8
                do
                    for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
                    do
                        mpirun -np $np $exe --number-of-timesteps           1025                \
                                            --number-of-time-grids          $nl                 \
                                            --coarsening-factors            $ml0 $ml1           \
                                            $lambdari                                           \
                                            --bound                         $bound              \
                                            --bound-on-level                1                   \
                                            --F-cycle                                           \
                                            --relaxation-scheme             F_relaxation        \
                                            --output-file                   TEST-$bound         &> /dev/null
                    done
                    val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
                    val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
                    absdiff=`python -c "print abs(($val0-$val1)/$val0)<$tol"`
                    numTests=$((numTests+1))
                    if [ "$absdiff" = "True" ]; then
                        numTestsPassed=$((numTestsPassed+1))
                    else
                        echo ">>>INFO: $nl-level F-cycle with F-relaxation. Test failed for parameters: m_0 = $ml0, m_1 = $ml1. $val0 != $val1"
                        numTestsFailed=$((numTestsFailed+1))
                    fi
                    if [[ $full != true ]]; then
                        break
                    fi
                done
                if [[ $full != true ]]; then
                    break
                fi
            done
        done
        if [[ $full != true ]]; then
            break
        fi
        # second, test odd coarsening factors
        testlambdar=$testdir/phi_evals_nt6561_nl6_ml3_nn11x11_L_stable_SDIRK$p\_c10_real.txt
        testlambdai=$testdir/phi_evals_nt6561_nl6_ml3_nn11x11_L_stable_SDIRK$p\_c10_imag.txt
        for lambdari in "--file-phi-real-eigenvalues $testlambdar" "--file-phi-complex-eigenvalues $testlambdar $testlambdai"
        do
            for ml0 in 3 9 27
            do
                for ml1 in 3 9
                do
                    for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
                    do
                        mpirun -np $np $exe --number-of-timesteps           730                 \
                                            --number-of-time-grids          $nl                 \
                                            --coarsening-factors            $ml0 $ml1           \
                                            $lambdari                                           \
                                            --bound                         $bound              \
                                            --bound-on-level                1                   \
                                            --F-cycle                                           \
                                            --relaxation-scheme             F_relaxation        \
                                            --output-file                   TEST-$bound         &> /dev/null
                    done
                    val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
                    val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
                    absdiff=`python -c "print abs(($val0-$val1)/$val0)<$tol"`
                    numTests=$((numTests+1))
                    if [ "$absdiff" = "True" ]; then
                        numTestsPassed=$((numTestsPassed+1))
                    else
                        echo ">>>INFO: $nl-level F-cycle with F-relaxation. Test failed for parameters: m_0 = $ml0, m_1 = $ml1. $val0 != $val1"
                        numTestsFailed=$((numTestsFailed+1))
                    fi
                done
            done
        done
    done
fi

if [[ $vcyc == true ]]; then
    echo ">>>INFO: Running 2-level cycle tests."
    ########################################################################################################################
    # compare 2-level cycle with F-/FCF-relaxation: inequality bound with inequality bound expression
    nl=2
    ########################################################################################################################
    # note: we are using arbitrary lambda's.
    #       thus, we cannot expect any reasonable values
    #       and in particular not: convergence estimate smaller than 1
    # for time integration orders
    for (( i = 0; i < ${#ps[@]}; ++i ))
    do
        p=${ps[i]}
        # for numbers of levels
        for r in F_relaxation FCF_relaxation
        do
            # first, test even coarsening factors
            testlambdar=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_real.txt
            testlambdai=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_imag.txt
            for lambdari in "--file-phi-real-eigenvalues $testlambdar" "--file-phi-complex-eigenvalues $testlambdar $testlambdai"
            do
                for ml0 in 2 4 8 16 32
                do
                    for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
                    do
                        mpirun -np $np $exe --number-of-timesteps           1025                \
                                            --number-of-time-grids          $nl                 \
                                            --coarsening-factors            $ml0                \
                                            $testlambdari                                       \
                                            --bound                         $bound              \
                                            --bound-on-level                1                   \
                                            --relaxation-scheme             $r                  \
                                            --output-file                   TEST-$bound         &> /dev/null
                    done
                    val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
                    val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
                    absdiff=`python -c "print abs(($val0-$val1)/$val0)<$tol"`
                    numTests=$((numTests+1))
                    if [ "$absdiff" = "True" ]; then
                        numTestsPassed=$((numTestsPassed+1))
                    else
                        echo ">>>INFO: $nl-level cycle with $r test failed for parameters: m_0 = $ml0. $val0 != $val1"
                        numTestsFailed=$((numTestsFailed+1))
                    fi
                    if [[ $full != true ]]; then
                        break
                    fi
                done
            done
            if [[ $full != true ]]; then
                continue
            fi
            # second, test odd coarsening factors
            testlambdar=$testdir/phi_evals_nt6561_nl6_ml3_nn11x11_L_stable_SDIRK$p\_c10_real.txt
            testlambdai=$testdir/phi_evals_nt6561_nl6_ml3_nn11x11_L_stable_SDIRK$p\_c10_imag.txt
            for lambdari in "--file-phi-real-eigenvalues $testlambdar" "--file-phi-complex-eigenvalues $testlambdar $testlambdai"
            do
                for ml0 in 3 9 27 81
                do
                    for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
                    do
                        mpirun -np $np $exe --number-of-timesteps           730                 \
                                            --number-of-time-grids          $nl                 \
                                            --coarsening-factors            $ml0                \
                                            $testlambdari                                       \
                                            --bound                         $bound              \
                                            --bound-on-level                1                   \
                                            --relaxation-scheme             $r                  \
                                            --output-file                   TEST-$bound         &> /dev/null
                    done
                    val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
                    val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
                    absdiff=`python -c "print abs(($val0-$val1)/$val0)<$tol"`
                    numTests=$((numTests+1))
                    if [ "$absdiff" = "True" ]; then
                        numTestsPassed=$((numTestsPassed+1))
                    else
                        echo ">>>INFO: $nl-level cycle with $r test failed for parameters: m_0 = $ml0. $val0 != $val1"
                        numTestsFailed=$((numTestsFailed+1))
                    fi
                done
            done
        done
    done
    echo ">>>INFO: Running 3-level V-cycle tests."
    ########################################################################################################################
    # compare 3-level V-cycle with F-/FCF-relaxation: inequality bound with inequality bound expression
    nl=3
    ########################################################################################################################
    # note: we are using arbitrary lambda's.
    #       thus, we cannot expect any reasonable values
    #       and in particular not: convergence estimate smaller than 1
    # for time integration orders
    for (( i = 0; i < ${#ps[@]}; ++i ))
    do
        p=${ps[i]}
        # for numbers of levels
        for r in F_relaxation FCF_relaxation
        do
            # first, test even coarsening factors
            testlambdar=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_real.txt
            testlambdai=$testdir/phi_evals_nt4096_nl6_ml2_nn11x11_L_stable_SDIRK$p\_c10_imag.txt
            for lambdari in "--file-phi-real-eigenvalues $testlambdar" "--file-phi-complex-eigenvalues $testlambdar $testlambdai"
            do
                for ml0 in 2 4 8 16
                do
                    for bound in error_l2_sqrt_expression_upper_bound error_l2_sqrt_upper_bound
                    do
                        mpirun -np $np $exe --number-of-timesteps           1025                \
                                            --number-of-time-grids          $nl                 \
                                            --coarsening-factors            $ml0 2              \
                                            $testlambdari                                       \
                                            --bound                         $bound              \
                                            --bound-on-level                1                   \
                                            --relaxation-scheme             $r                  \
                                            --output-file                   TEST-$bound         &> /dev/null
                    done
                    val0=`cat max_TEST-error_l2_sqrt_expression_upper_bound.txt`
                    val1=`cat max_TEST-error_l2_sqrt_upper_bound.txt`
                    absdiff=`python -c "print abs(($val0-$val1)/$val0)<$tol"`
                    numTests=$((numTests+1))
                    if [ "$absdiff" = "True" ]; then
                        numTestsPassed=$((numTestsPassed+1))
                    else
                        echo ">>>INFO: $nl-level V-cycle with $r test failed for parameters: m_0 = $ml0, m_1 = $ml1. $val0 != $val1"
                        numTestsFailed=$((numTestsFailed+1))
                    fi
                    if [[ $full != true ]]; then
                        break
                    fi
                done
            done
            if [[ $full != true ]]; then
                continue
            fi
            # second, test odd coarsening factors
            # todo
        done
    done
fi

end=`date +%s`

runtime=$((end-start))

if [[ $numTests == 0 ]]; then
    echo ">>>INFO: No tests were run. Please check arguments."
else
    echo ""
    echo ">>>INFO: Passed tests: $numTestsPassed / $numTests"
    echo ""
    echo ">>>INFO: Elapsed time: $runtime s"
    echo ""
fi

cd $here
