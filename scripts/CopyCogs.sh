#!/bin/bash

OLDIFS=$IFS
IFS=,
while read cluster taxaid 
do
    echo $cluster
#    echo $taxaid

    mkdir "Strain_$taxaid/ClusterCogs"
    cp Simulation/Split/${cluster}/${cluster}_core.cogs Strain_$taxaid/ClusterCogs
done < ClusterSpecies.txt
IFS=$OLDIFS
