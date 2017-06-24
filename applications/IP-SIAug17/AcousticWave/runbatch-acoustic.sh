#!/bin/bash

parameps=1e-3
paramks='5e-6 2e-6 9e-7 5e-7 2e-7 9e-8 5e-8'

echo 'eps='$parameps

for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 15 python exple_acousticinversion.py $paramk $parameps > 'exple_acousticinversion/k'$paramk'_e'$parameps
done

echo Bash completed
