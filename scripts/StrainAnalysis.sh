#!/bin/bash

speciesId=$1

./MoveStrains.pl $speciesId < assembly_summary.txt

cd Strain_${speciesId}

echo "Run prodigal"

for file in GCF*/*fna
do
    echo $stub

    stub=${file%.fna}
    
    prodigal -p meta -i $file -a ${stub}.faa -d ${stub}.fas -f gff -o ${stub}.gff > $stub.out&

done

wait

echo "Perform RPS blast"

for file in GCF*/*faa
do
    stub=${file%.faa}
    rpsblast -outfmt "6 qseqid sseqid evalue pident score qstart qend sstart send length slen" -max_target_seqs 1 -evalue 0.00001 -db /mnt/data-chris/chris/Databases/rpsblast_db/Cog -num_threads 16 -query $file > ${stub}.out&
done

wait

echo "Assign COGs"

for dir in GCF*/
do
    cd $dir
    echo $dir
    for file in *.out
    do
        stub=${file%.out}
        python ~/bin/ExtractCogsNative.py -b $file --cdd_cog_file ~/Installed/CONCOCT/scgs/cdd_to_cog.tsv > ${stub}.cogs&
    done

    cd ..
done

wait

~/bin/StrainMap.pl $speciesId < ../assembly_summary.txt > strain_map.csv

~/bin/ExtractCOGs.pl ~/gpfs/Databases/scg_cogs_min0.97_max1.03_unique_genera_LN.txt strain_map.csv

mkdir SCGs

mv *ffn SCGs

cat SCGs/*ffn > temp.fa
grep ">" temp.fa | sed 's/>//g' | sed 's/_COG.*//' | sort | uniq > Hap.txt
~/bin/CombineH.pl Hap.txt SCGs/*ffn > SCGs.fa
mafft SCGs.fa > SCGs.gfa

FastTreeMP -nt -gtr < SCGs.gfa > SCGs.tree

python ~/bin/IdentityH.py -i SCGs.gfa > IdentH.txt

mkdir GeneClusters

cd GeneClusters

cat ../GCF_00*/*fas > All.fas

vsearch --cluster_fast All.fas --id 0.95 --centroids Cen.fa

cd ..

for file in GCF_*/*fna
do
    stub=${file%.fna}
    echo $stub
    vsearch -usearch_global ${stub}.fas -db GeneClusters/Cen.fa -strand plus -id 0.95 -uc ${stub}_fas_map.uc
done

~/bin/AddOverlap.pl strain_map.csv < IdentH.txt > IdentHG.csv

cd ..
