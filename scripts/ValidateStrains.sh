#!/bin/bash

: "${DESMAN? Need to set DESMAN}"

while IFS=, read cluster strain tp total predicted recall precision
do

    if (( $predicted > 0 )); then
        cd ${cluster}_scg
        printf "$cluster,$strain,"
        python $DESMAN/scripts/resolvenhap.py -d 0.05 -m 0.10 -f 0.05 ${cluster}_scg 
        cd ..
    fi
done < VarResults_V.csv

 
