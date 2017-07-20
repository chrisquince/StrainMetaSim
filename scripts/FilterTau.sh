#!/bin/bash
: "${METASIMPATH?Need to set METASIMPATH}"
OLDIFS=$IFS
IFS=,
while read cluster taxaid 
do
    echo $cluster
#    echo $taxaid

    if [ -f  Strain_${taxaid}/${cluster}_core_tau_map.csv ]; then
    $METASIMPATH/FilterTau.pl Simulation/SCG_Analysis/${cluster}_scg/${cluster}_scgcogf.csv < Strain_${taxaid}/${cluster}_core_tau_map.csv > Strain_${taxaid}/${cluster}_core_tau_mapF.csv
    fi
done < ClusterSpecies.txt
IFS=$OLDIFS
