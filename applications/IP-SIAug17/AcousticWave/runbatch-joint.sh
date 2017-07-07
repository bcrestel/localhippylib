#!/bin/bash

paramks='0.0 1e-12 1e-10 1e-8 1e-6 1e-4'


for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 8 python joint_poisson-acoustic.py $paramk > 'joint_poisson-acoustic/joint4Hz_TV+CG_k'$paramk'_e1e-3.out' &
done

echo Bash completed for joint
