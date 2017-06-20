#!/bin/bash

parameps=1e-1
paramks='1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8'

echo 'eps='$parameps

for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 15 python exple_acousticinversion.py $paramk $parameps > 'exple_acousticinversion/k'$paramk'_e'$parameps
done

echo Bash completed
