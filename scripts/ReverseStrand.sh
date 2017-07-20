#!/bin/bash

OLDIFS=$IFS
IFS=,
while read species nstrains 
do
    echo $species
#    echo $taxaid

    if (( $nstrains > 1 )); then
        cd Strain_${species}
        for cogFile in ClusterCogs/*_core.cogs 
        do 
            if [ -e $cogFile ]; then
                basename=${cogFile##*/} 
                cluster=${basename%_core.cogs}
                cp -r Select_SCGs Select_SCGs_${cluster}
                ../ReverseStrand.pl $cogFile Select_SCGs_${cluster}
                ../TauFastaC.pl Select_SCGs_${cluster}
 
                ../CombineTau.pl Select_SCGs_${cluster} > ${cluster}_core_tau.csv
                ../MapCogBack.pl $cogFile < ${cluster}_core_tau.csv > ${cluster}_core_tau_map.csv 
              echo $cluster
            fi
        done
        cd .. 
    fi

done < StrainCount.csv
IFS=$OLDIFS
