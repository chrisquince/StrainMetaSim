#!/bin/bash
: "${METASIMPATH?Need to set METASIMPATH}"
for file in */*_core_tau_mapF.csv
do
    stub=${file%.csv}
    stub2=${file%\/Cluster*tau_mapF.csv}
    printf "$stub2,"
    python $METASIMPATH/scripts/Degenerate.py $file ${stub}U.csv

done
