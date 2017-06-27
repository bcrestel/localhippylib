#!/bin/bash

parameps=1e-3
paramks='5e-7 4e-7 3e-7 2e-7 1e-7 9e-8 8e-8'

echo 'eps='$parameps

for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 36 python exple_acousticinversion.py $paramk $parameps > 'exple_acousticinversion/k'$paramk'_e'$parameps
done

echo Bash completed for acoustic
