#!/bin/bash

paramks='6e-8 7e-8 8e-8 9e-8 1e-7 2e-7 4e-7'


for paramk in $paramks
do
    echo 'k='$paramk
    nohup mpirun -n 8 python joint_poisson-acoustic.py $paramk > 'joint_poisson-acoustic/joint4Hz_NN_k'$paramk'_e1e-3.out' &
done

echo Bash completed for joint
