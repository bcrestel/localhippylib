#!/bin/bash

parameps=1e-3
paramks='5e-9 7e-9 1e-8 2e-8 3e-8 5e-8'

for paramk in $paramks
do
    echo 'paramk='$paramk
    nohup python exple_poisson.py $paramk $parameps > 'exple_poisson/k'$paramk'_e'$parameps'.out'
done

echo Bash completed for poisson
