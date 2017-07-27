#!/bin/bash
#_core_tau_map.csv
#Cluster132,Strain_2095,17,18,18,0.944444444444,0.944444444444

while IFS=, read cluster strain tp total predicted recall precision
do

    if (( $predicted > 0 )); then
        cd ${cluster}_scg
        printf "$cluster,$strain,"
        python $DESMAN/scripts/resolvenhap.py -d 0.05 -m 0.10 -f 0.05 ${cluster}_scg 
        cd ..
    fi
    #../DESMAN/scripts/validateSNP2.py $cluster/$cluster_scg ../../$strain/$cluster_core_tau_map.csv

done < VarResults_V.csv

 
