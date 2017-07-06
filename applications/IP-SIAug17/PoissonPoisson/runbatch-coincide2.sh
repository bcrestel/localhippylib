#!/bin/bash

paramks='7e-7 5e-7 1e-7 1e-6 5e-6'

for paramk in $paramks
do
    echo 'paramk='$paramk
    mpirun -n 16 python joint-reconstructions-coincide2.py $paramk > 'output/TVPD+'$paramk'NCG_eps1e-6.out'
done

echo Bash completed for coincide2
