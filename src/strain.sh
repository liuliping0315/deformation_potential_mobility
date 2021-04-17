#!/bin/bash
for j in 1 2
do      
    for k in 0.990 0.995 1.000 1.005 1.010
    do
        echo start $j $k
        cp -ra 0 ${j}_${k}
        cd ${j}_${k}
        strain.py ${j} ${k}
        \mv poscar_${k} POSCAR
        cd ..
    done
done

