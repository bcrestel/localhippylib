#!/bin/bash

parameps=1e-3
paramks='1e-6 5e-7 3e-7 2e-7 1e-7 9e-8 8e-8 5e-8 1e-8'

for paramk in $paramks
do
    echo 'paramk='$paramk
    nohup python exple_poisson.py $paramk $parameps > 'exple_poisson/k'$paramk'_e'$parameps &
done

echo Bash completed for poisson
