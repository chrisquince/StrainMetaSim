#!/bin/bash
#_core_tau_map.csv

: "${METASIMPATH?Need to set METASIMPATH}"

while IFS=, read cluster strain
do
    tauFile=../../Strain_${strain}/${cluster}_core_tau_mapFU.csv
    if [ -f $tauFile ]; then
        printf "$cluster,$strain,"
        python $METASIMPATH/scripts/var_validate.py ${cluster}_scg/${cluster}_scgsel_var.csv $tauFile 
    else

        a=($(wc ${cluster}_scg/${cluster}_scgsel_var.csv))
        let nvar=${a[0]}-1
        echo "NV,$cluster,$strain,$nvar"

    fi
        #    echo "I got:$cluster|$strain"
done < ClusterSpecies.txt

 
