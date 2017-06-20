#!/bin/bash

parameps=1e-4
paramks='5e-3 4e-3 3e-3 2e-3 1e-3 9e-4 8e-4 7e-4 6e-4 5e-4'

for paramk in $paramks
do
    echo 'paramk='$paramk
    nohup python exple_poisson.py $paramk $parameps > 'exple_poisson/k'$paramk'_e'$parameps &
done

echo Bash completed
