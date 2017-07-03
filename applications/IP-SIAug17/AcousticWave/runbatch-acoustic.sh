#!/bin/bash

parameps=1e-3
paramks='9e-7 8e-7 7e-7 6e-7'

echo 'eps='$parameps

for paramk in $paramks
do
    echo 'k='$paramk
    mpirun -n 30 python exple_acousticinversion.py $paramk $parameps > 'exple_acousticinversion/k'$paramk'_e'$parameps
done

echo Bash completed for acoustic
