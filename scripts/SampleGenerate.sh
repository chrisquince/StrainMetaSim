#!/bin/bash

count=0
for ((n=0; n<96; n++))
{
     python ./SampleGenerate.py Select_config.json Simulation Simulation/coverage.tsv -n Reads -s $n > log_$n.out& 

    if [ $count -eq 48 ] ; then

        wait
        count=0
    fi

    let count=count+1
}


