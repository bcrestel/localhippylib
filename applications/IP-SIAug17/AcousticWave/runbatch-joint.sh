#!/bin/bash

paramks='2e-8 1e-8 3e-8 5e-8 5e-9'


for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 8 python joint_poisson-acoustic.py $paramk > 'joint_poisson-acoustic/joint4Hz_NN_k'$paramk'_e1e-3.out'
done

echo Bash completed for joint
