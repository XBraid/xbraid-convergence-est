#!/bin/sh

exe=./main
np=4

Nl=1025
ml=2
#cf0=2
cf=2
#relax=F_relaxation
bound=error_l2_tight_twogrid_upper_bound
boundOnLevel=0
#rk=A_stable_LobattoIIIA_order2
rmin=-100.0
rmax=100.0
imin=-100.0
imax=100.0
rnum=1051
inum=1051

for cf0 in 2 4 8 16
do
    for relax in F_relaxation FCF_relaxation
    do
        for rk in A_stable_LobattoIIIA_order2 A_stable_LobattoIIIA_order4 A_stable_LobattoIIIA_order6 A_stable_LobattoIIIB_order2 A_stable_LobattoIIIB_order4 A_stable_LobattoIIIB_order6 A_stable_LobattoIIIB_order8 A_stable_Gauss_order2 A_stable_Gauss_order4 A_stable_Gauss_order6 A_stable_SDIRK2 A_stable_SDIRK3 A_stable_SDIRK4 L_stable_SDIRK1 L_stable_SDIRK2 L_stable_SDIRK3 L_stable_SDIRK4 L_stable_SDIRK5 L_stable_RadauIIA_order1 L_stable_RadauIIA_order3 L_stable_RadauIIA_order5 L_stable_LobattoIIIC_order2 L_stable_LobattoIIIC_order4 L_stable_LobattoIIIC_order6 L_stable_LobattoIIIC_order8 LobattoIIICast_order2 LobattoIIICast_order4 LobattoIIICast_order6 LobattoIIICast_order8
        do
            outname=range_$rmin\_$rmax\_$imin\_$imax\_num_$rnum\_$inum\_ml_$ml\_$relax\_cf0_$cf0\_bound_$bound\_level_$bound\_rk_$rk
            echo "$outname"
            mpirun -np $np $exe --number-of-timesteps $Nl --number-of-time-grids $ml --coarsening-factors $cf0 \
                                --sample-complex-plane $rmin $rmax $imin $imax --complex-plane-sample-size $rnum $inum \
                                --bound $bound --bound-on-level $boundOnLevel \
                                --relaxation-scheme $relax \
                                --runge-kutta-method $rk \
                                --output-file $outname
            sleep 1
            mv dteta_l0.txt dteta_l0_$outname.txt
            sleep 1
            echo " ==== rendering results"
            python3 utils/heatmap.py dteta_l0_$outname.txt $outname.txt NAN NAN NAN NAN 0 1
            sleep 1
            mv tst.png range_$rmin\_$rmax\_$imin\_$imax\_num_$rnum\_$inum\_ml_$ml\_$relax\_cf0_$cf0\_bound_$bound\_level_$bound\_rk_$rk.png
            echo " ==== done"
        done
    done
done