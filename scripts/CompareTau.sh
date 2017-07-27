#!/bin/bash

#Cluster132,4,2,3,0.0,Cluster132_scg_4_3/Filtered_Tau_star.csv

: "${METASIMPATH?Need to set METASIMPATH}"

: "${DESMAN?Need to set DESMAN}"

while IFS=, read cluster strain G H R acc prediction
do
    stub=${prediction%_star.csv}

    printf "$cluster,$strain,"

    a=($(wc ${cluster}_scg/${cluster}_scgsel_var.csv))
    let nvar=${a[0]}-1

    printf "$nvar,"

    python $DESMAN/scripts/validateSNP3.py ${cluster}_scg/${stub}_starR.csv $METASIMPATH/ComplexStrainSim/Strains/Strain_${strain}/${cluster}_core_tau_mapFU.csv

done < ValidateStrains.csv

 
