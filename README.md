# StrainMetaSim

This repository contains scripts and downloads to data necessary to recreate the `complex strain' in silico mock community from Quince et al "DESMAN: a new tool to decipher strain-level divergence de novo from metagenomes".

The strategy used here (and some of the code) was heavily based on the excellent and far more general [meta-sweeper](https://github.com/cerebis/meta-sweeper) 

1. Lets create a directory for working in and set PATH to scripts:
    ```
        mkdir ComplexStrainSim
        cd ComplexStrainSim
        export METASIMPATH=/MyPath/StrainMetaSim/
    ```
 

2. Then download the processed data for 100 species with multiple strains obtained from the NCBI.
    ```
        wget https://complexstrainsim.s3.climb.ac.uk/Strains.tar.gz
    ```
These directories were created from NCBI genomes using a script StrainAnalysis.sh which is included in the repo. If you want to prepare your own species for inclusion in the simulation then you should run this script as follows:
    ```
        StrainAnalysis.sh taxa_id
    ```
This will require that you are in directory with NCBI genomes and the file assembly_summary.txt best to contact me if you want to do this.

3. Extract out species directories:
    ```
        tar -xvzf Strains.tar.gz
        cd Strains
    ```

4. You will now have 100 directories one for each species, each of which contains detailed information on the strains available from the NCBI for that strain. 

5. The actual simulation is configured through a json file. Copy that into the Strains directory:
    ```
        cp ../../config/Select_config.json .
    ```
    and view it in an editor. The parameters are fairly self explanatory:
    ```
        "reads": {
            "no_Reads": 12500000,

            "length": 150,

            "insert_length": 300,

            "insert_sd": 10
        },
    ```
    Configures the Illumina paired end reads, the number for each sample, their length, and insert characteristics.
    ```  
        "no_Samples" : 96,
    ```
    This is the number of samples to generate. 
    ```
        "species_dist_Params": {
            "mean_log_Mean": 1.0,

            "sd_log_Mean": 0.25,

            "k_log_Sd": 1.0,

            "theta_log_Sd": 1.0,

            "beta": 1.0,

            "alpha": 1.0
        },
    ```
    These are the parameters for the Normal (mean_log_Mean, sd_log_Mean) and Gamma (k_log_Sd, theta_log_Sd) 
    distributions from which we obtain the log-normal parameters (mean log, sd log) for each species. Beta is 
    an optional parameter to remove some species from samples entirely and alpha is the parameter for the symmetric 
    Dirichlet used to determine the relative strain frequencies within a species. 

    ```
        "no_Species" : 100,

        "species": [
        {
            "dir": "Strain_35814/",
            "nStrains": "4"
        },
    ```
    Then we define the total number of species in the simulation, give their directory names and the number of strains to generate 
    for them.

6. Now we begin by generating the coverages required for each strain:
    ```
        python $METASIMPATH/scripts/CoverageGenerate.py Select_config.json Simulation
    ```
    Where Simulation is the output directory name which will be created if not present.         This will generate one fasta file named foo.tmp for each strain in the collection, this contains the sequences for that strain. The strains that have been selected for the simulation are detailed in the tab separated *select.tsv* file:
    ```
    head -n 5 select.tsv 
    35814	1247648_0	4	GCF_000598125.1
    35814	1172205_1	1	GCF_000765395.1
    35814	1247646_0	4	GCF_000572015.1
    35814	1266729_0	2	GCF_000341465.1
    285	399795_0	1	GCF_000168855.1
    ```
    This simply gives taxaid strainid no_seqs accession. Then the file coverage.tsv gives the coverage of each strain in each sample:
    ```
    head -n 5 coverage.tsv
    0	663	1219076_0	0.212580070954	0.000583505311272
    0	2095	1124992_0	0.55214734603	0.000299582781855
    0	2095	40480_2	0.259023748867	0.000140534958164
    0	2096	710127_0	3.80894852936e-06	2.05744163762e-09
    0	2096	1159198_0	3.50903469472e-07	1.78537626634e-10
    ```
    The first column in the sample index, followed by taxaid, strainid, coverage, and relative frequency.
    
7. Now we will generate the actual reads basically by running the [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/) read simulator for each strain in each sample separately. We do this through the script
*SampleGenerate.py* which we run as follows:
```
python ./SampleGenerate.py Select_config.json Simulation Simulation/coverage.tsv -n Reads -s 0
```
where -s gives the sample index to be generated. This is trivially parallelisable so we run all the samples simultaneously with the following shell script:
```
    cp $METASIMPATH/scripts/SampleGenerate.sh .
    cp $METASIMPATH/scripts/SampleGenerate.py .
    ./SampleGenerate.sh
```
This will generate 96 paired end samples named Reads.0.r1.fq.gz, Reads.0.r2.fq.gz,...,Reads.95.r1.fq.gz, Reads.95.r2.fq.gz

We will move these into a separate reads directory along with genomes:
```
mkdir Reads
mv Reads*gz Reads
mkdir Genomes
mv *tmp Genomes
```

### Downloading reads (alternative step to above)

If you cannot get the above simulation to work or just want to go ahead with the DESMAN 
example then all the reads are downloadable from the URLs contained in this file [read urls](Results/complexmock.txt). Create a directory structure as above if you have not already i.e. starting from the repo dir: 
```
cd  $METASIMPATH
mkdir -p ComplexStrainSim/Strains/Simulation
cd ComplexStrainSim/Strains/Simulation
```

Now download reads:
```
mkdir Reads
cd Reads
while read line
do
    wget $line
done < $METASIMPATH/Results/complexmock.txt
cd ..
```

You will also need the genomes and parameter files for the simulation:
```
wget https://complexmock.s3.climb.ac.uk/Genomes.tar.gz
tar -xvzf Genomes.tar.gz
wget https://complexmock.s3.climb.ac.uk/select.tsv
wget https://complexmock.s3.climb.ac.uk/coverage.tsv
```

You can then proceed with the rest of the analysis below.

## Running DESMAN on the complex mock

We now describe how to perform a complete analysis binning and resolving strains on 
this synthetic community. Some of these steps are time consuming so we provide the option to download the output instead. This assumes that DESMAN and CONCOCT are installed and their paths 
set to the variables DESMAN and CONCOCT respectively e.g. (changing paths to your 
system set-up):

```
export CONCOCT=/mnt/gpfs/chris/repos/CONCOCT/
export DESMAN=/mnt/gpfs/chris/repos/DESMAN/
```

We will also create a new variable pointing to our current working dir for all this analysis:
```
export METASIMPATHWD=$METASIMPATH/ComplexStrainSim/Strains/Simulation
```

The first step in the analysis is to assemble the reads. 

### Assembly

We assembled the reads using MEGAHIT 1.1.1 and default parameters:
```
ls Reads/*r1*gz | tr "\n" "," | sed 's/,$//' > r1.csv
ls Reads/*r2*gz | tr "\n" "," | sed 's/,$//' > r2.csv
megahit -1 $(<r1.csv) -2 $(<r2.csv) -t 96 -o Assembly > megahit1.out
```

Then cut up contigs and index for BWA:

```
cd Assembly
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m final.contigs.fa > final_contigs_c10K.fa
bwa index final_contigs_c10K.fa
cd ..
```

If you do not have time to run the assembly the results can be downloaded here:
```
wget https://complexmockresults.s3.climb.ac.uk/Assembly.tar.gz
tar -xvzf Assembly.tar.gz
```

### Mapping

Then we map reads onto the contig fragments using BWA-mem:
```
mkdir Map

for file in Reads/*.r1.fq.gz
do

    stub=${file%.r1.fq.gz}

    base=${stub##*/}

    echo $base

    bwa mem -t 96 Assembly/final_contigs_c10K.fa $file ${stub}.r2.fq.gz > Map/$base.sam
done
```

And calculate coverages:
```
python ~/bin/Lengths.py -i Assembly/final_contigs_c10K.fa | tr " " "\t" > Assembly/Lengths.txt

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
    (samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam; samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam; bedtools genomecov - ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt)&
done
```

Collate coverages together:

```
for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv&
done

$DESMAN/scripts/Collate.pl Map > Coverage.csv
```

### Binning

We are going to use a more efficient development branch of CONCOCT. This can be checked out and installed as follows:

```
    git clone git@github.com:BinPro/CONCOCT.git
    cd CONCOCT
    git fetch
    git checkout SpeedUp_Mp
    sudo python ./setup.py install
```

Now we can run CONCOCT:
```

    mkdir Concoct

    mv Coverage.csv Concoct

    cd Concoct

    tr "," "\t" < Coverage.csv > Coverage.tsv

    concoct --coverage_file Coverage.tsv --composition_file ../Assembly/final_contigs_c10K.fa -t 96 > concoct.out

```

Find genes using prodigal:
```
    cd ..
    
    mkdir Annotate

    cd Annotate/

    python $DESMAN/scripts/LengthFilter.py ../Assembly/final_contigs_c10K.fa -m 1000 >     final_contigs_gt1000_c10K.fa

    prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d     final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff > p.out
```

Assign COGs change the -c flag which sets number of parallel processes appropriately:
```
    export COGSDB_DIR=/mydatabase_path/rpsblast_db
    $CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 32 -r 1
```

We are also going to refine the output using single-core gene frequencies. First we calculate scg frequencies on the CONCOCT clusters:
```
cd ../Concoct
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_gt1000.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv
```

Then we need to manipulate the output file formats slightly:
```
sed '1d' clustering_gt1000.csv > clustering_gt1000_R.csv
cut -f1,3- < clustering_gt1000_scg.tsv | tr "\t" "," > clustering_gt1000_scg.csv
$CONCOCT/scripts/Sort.pl < clustering_gt1000_scg.csv > clustering_gt1000_scg_sort.csv
```

Then we can run the refinement step of CONCOCT:
```
concoct_refine clustering_gt1000_R.csv original_data_gt1000.csv clustering_gt1000_scg_sort.csv > concoct_ref.out
```

This should result in 70-75 clusters with 75% single copy copy SCGs:
```
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_refine.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_refine_scg.tsv
```

If the above fails or time prevents you from then the CONCOCT results can be downloaded directly:

```
wget https://complexmockresults.s3.climb.ac.uk/Concoct_Res.tar.gz
tar -xvzf Concoct_Res.tar.gz 
mv Concoct_Res Concoct
```


### Get nucleotide frequencies on target bins

For this part of the analysis we will require an additional github repository of scripts
[MAGAnalysis](https://github.com/chrisquince/MAGAnalysis). Install as follows but into 
you own repositories directory:
```
cd /mnt/gpfs/chris/repos/
git clone https://github.com/chrisquince/MAGAnalysis
cd $COMPLEXSIMWD
export MAGANALYSIS=/mnt/gpfs/chris/repos/MAGAnalysis
```  

Now we assign COGs and genes to contigs using one of these scripts:
```
cd $COMPLEXSIMWD/Annotate
python $DESMAN/scripts/ExtractCogs.py -b final_contigs_gt1000_c10K.out --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.cogs
python /mnt/gpfs/chris/repos/DESMAN/scripts/ExtractGenes.py -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.genes
cd ..
```

Return to the analysis directory and create a new directory to bin the contigs into:
```
cd $COMPLEXSIMWD
mkdir Split
$DESMAN/scripts/SplitClusters.pl ../Annotate/final_contigs_gt1000_c10K.fa ../Concoct/clustering_refine.csv
$METASIMPATH/scripts/SplitCOGs.pl ../Annotate/final_contigs_gt1000_c10K.cogs ../Concoct/clustering_refine.csv
$METASIMPATH/scripts/SplitGenes.pl ../Annotate/final_contigs_gt1000_c10K.genes ../Concoct/clustering_refine.csv
cd ..
```

Now we want to select those clusters that have 75% of SCGs in single copy using R:
```
cd Concoct
R
```

Then run these R commands:
```
> scg_ref <- read.table("clustering_refine_scg.tsv",header=TRUE,row.names=1)
> scg_ref <- scg_ref[,-1]
> scg_ref <- scg_ref[,-1]
> scg_75 <- scg_ref[rowSums(scg_ref==1)/36 > 0.75,]
> write.csv(scg_75,"scg_75.csv",quote=FALSE)
> q()
```

Back in bash:
```
cut -f1 -d"," scg_75.csv | sed 's/^/Cluster/' | sed '1d' > Cluster75.txt
cd ..
```

Now we can split up bam files by each cluster in turn:
```
mkdir SplitBam

while read -r cluster 
do
    grep ">" Split/${cluster}/${cluster}.fa | sed 's/>//g' > Split/${cluster}/${cluster}_contigs.txt
    $METASIMPATH/scripts/AddLengths.pl Annotate/final_contigs_gt1000_c10K.len < Split/${cluster}/${cluster}_contigs.txt > Split/${cluster}/${cluster}_contigs.tsv
    mkdir SplitBam/${cluster}

    for bamfile in Map/*.mapped.sorted.bam
    do
        stub=${bamfile#Map\/}
        stub=${stub%.mapped.sorted.bam}
        
        samtools view -bhL Split/${cluster}/${cluster}_contigs.tsv $bamfile > SplitBam/${cluster}/${stub}_Filter.bam&

    done 
    wait    
done < Concoct/Cluster75.txt 
```

and run bam-readcount:
```
while read line
do
    mag=$line

    echo $mag

    cd SplitBam
    cd ${mag}

    cp ../../Split/${mag}/${mag}_contigs.tsv ${mag}_contigs.tsv
    samtools faidx ../../Split/${mag}/${mag}.fa

    echo "${mag}_contigs.tsv"
    mkdir ReadcountFilter
    for bamfile in *_Filter.bam
    do
        stub=${bamfile%_Filter.bam}
            
        echo $stub

        (samtools index $bamfile; bam-readcount -w 1 -q 20 -l "${mag}_contigs.tsv" -f ../../Split/${mag}/${mag}.fa $bamfile 2> ReadcountFilter/${stub}.err > ReadcountFilter/${stub}.cnt)&
    
     done
     cd ..
     cd ..
     wait
    
done < Concoct/Cluster75.txt

```

Then we want to select core cogs from each cluster

```
while read -r cluster 
do
    echo $cluster
    $DESMAN/scripts/SelectContigsPos.pl $METASIMPATH/config/cogs.txt < Split/${cluster}/${cluster}.cog > Split/${cluster}/${cluster}_core.cogs
done < Concoct/Cluster75.txt
```

Then we can get the base counts on these core cogs:

```
#!/bin/bash

mkdir Variants
while read -r cluster 
do
    echo $cluster
    cd ./SplitBam/${cluster}/ReadcountFilter
    gzip *cnt
    cd ../../..
    python $DESMAN/scripts/ExtractCountFreqGenes.py Split/${cluster}/${cluster}_core.cogs ./SplitBam/${cluster}/ReadcountFilter --output_file Variants/${cluster}_scg.freq > Variants/${cluster}log.txt&

done < Concoct/Cluster75.txt 
``` 

The results of all the above analysis the nucleotide frequencies on core genes in each cluster can be downloaded if the above fails and you can proceed to the variant detection:

```
wget https://complexmockresults.s3.climb.ac.uk/Variants.tar.gz
tar -xvzf Variants.tar.gz
```

### Variant detection

The first step of DESMAN is to find variant positions on the single-copy core genes. This is done using the Variant_Filter.py script. We run this for each cluster in parallel below adapt this to your machines capabilities:

```
cd $COMPLEXSIMWD

mkdir SCG_Analysis

for file in ./Variants/*.freq
do
    stub=${file%.freq}
    stub=${stub#./Variants\/}

    echo $stub
    mkdir SCG_Analysis/$stub
    
    cp $file SCG_Analysis/$stub
    cd SCG_Analysis/$stub    

    python $DESMAN/desman/Variant_Filter.py ${stub}.freq -p -o $stub -m 1.0 -f 25.0 -c -sf 0.80 -t 2.5 > ${stub}.vout&
    
    cd ../..
done
```

The parameters passed to Variant_Filter.py above are:
1. *-p* : Use 1d optimisation for minor variant frequency, slower but more sensitive 
2. *-o $stub* : Output file stub name 
3. *-m 1.0* : Minimum coverage for a sample not to be filtered
4. *-f 25.0* : Cut-off for initial variant filtering
5. *-c* : Filter genes/COGs by median coverage
6. *-sf 0.80* : Fraction of samples *i.e.* 80% that need to pass coverage filter for gene to be **not** filtered
7. *-t 2.5* : Coverage divergence from median for COG to be filtered

### Resolve haplotypes on core genes

We use the positions identified as variants above for the haplotype inference. The script below will run desman for each cluster in turn, for between 1 and 7 haplotypes for 10 different replicates. It then pauses, so this code will use 70 threads at a time, adjust as is appropriate for your computing environment:
```
#!/bin/bash

for dir in Cluster*/ 
do
    cd $dir
    echo $dir
    stub=${dir%\/}
    echo $stub    
    varFile=${stub}sel_var.csv

    eFile=${stub}tran_df.csv
    

    for g in 1 2 3 4 5 6 7 
    do
        for r in 0 1 2 3 4 5 6 7 8 9
            do
            echo $g
                (desman $varFile -e $eFile -o ${stub}_${g}_${r} -g $g -s $r -m 1.0 > ${stub}_${g}_${r}.out)&
            done
    done

    wait
    cd ..
done
```

If you fail to run the above variant detection and haplotype resolution then the 
results can be downloaded from here:

```
wget https://complexmockresults.s3.climb.ac.uk/SCG_Analysis.tar.gz
tar -xvzf SCG_Analysis.tar.gz
```

# Validation of variant and haplotype detection results

## Assign contigs to genomes

First we need to index the sorted per sample bam files:
```
cd $COMPLEXSIMWD
for file in Map/*mapped.sorted.bam
do     
    stub=${file%.bam}     
    stub2=${stub#Map\/}     
    echo $stub     
    samtools index $file 
done
```

Make a directory and concatenate the individual genome sequences into a single file:

```
mkdir AssignContigs
cd AssignContigs
cat ../Genomes/*tmp > AllGenomes.fasta
```

Then we run a slightly adjusted version of the contig_read_count_per_genome script from the DESMAN main repo:

```
python $COMPLEXSIM/scripts/contig_read_count_per_genomeM.py ../Assembly/final_contigs_c10K.fa AllGenomes.fasta ../Map/*mapped.sorted.bam > final_contigs_c10K_counts.tsv
```

This make take while so run in the background with screen etc. We then process this file to compute for each contig the species and strain the majority of reads derive from.

```
python ./MapCounts.py ../Genomes ../select.tsv final_contigs_c10K_counts.tsv
$COMPLEXSIM/scripts/Filter.pl < Species.csv > Contig_Species.csv
$COMPLEXSIM/scripts/Filter.pl < Strain.csv > Contig_Strains.csv 
```

Can now look at accuracy of CONCOCT clustering against species assignments:
```
$CONCOCT/scripts/Validate.pl --cfile=../Concoct/clustering_refine.csv --sfile=Contig_Species.csv --ffile=../Annotate/final_contigs_gt1000_c10K.fa
```

The result should look something like:
```
N	M	TL	S	K	Rec.	Prec.	NMI	Rand	AdjRand
74518	74485	4.0917e+08	100	144	0.848982	0.971654	0.956139	0.994957	0.812190
```

Now use R to map clusters onto species:
```
R
```
Run the following R comands:
```
>Conf <- read.csv("Conf.csv",header=TRUE,row.names=1)
>Conf <- t(Conf)
>scg <- read.table("../Concoct/clustering_refine_scg.tsv",header=TRUE,row.name=1)
>scg <- scg[,-1]
>scg <- scg[,-1]
>rownames(scg) <- gsub("^","D",rownames(scg))
>ConfR <- Conf[rownames(scg),]
>Conf75 <- ConfR[rowSums(scg==1)/36 > 0.75,]
>Conf75P <- Conf75/rowSums(Conf75)
>clust_species <- cbind.data.frame(Species=colnames(Conf75P)[apply(Conf75P,1,which.max)],F=apply(Conf75P, 1, max))
>write.csv(clust_species,"clust_species.csv",quote=FALSE)
>q()
```

The above processing can also be downloaded to save time:
```
wget https://complexmockresults.s3.climb.ac.uk/SCG_Analysis.tar.gz
tar -xvzf SCG_Analysis.tar.gz
```

We are now going to create for each species that maps onto a cluster what the expected variants should be from the genomes themselves. We move to the directory two up in the hierarchy *../StrainMetaSim/ComplexStrainSim/Strains*:

```
cd ../..
cut -f1 < Simulation/select.tsv | sort | uniq -c | sed 's/^[ \t]*//;s/[ \t]*$//' | awk '{ print $2 " " $1}' | tr " " "," > StrainCount.csv
```

and count number of strains in each species. Then we copy the mapping file and reformat it:

```
cp Simulation/AssignContigs/clust_species.csv .
cut -f1,2 -d"," clust_species.csv | sed 's/D/Cluster/' > ClusterSpecies.txt
```

Now we copy in the core cog files from each cluster to the correct species and use these we run a script to find the true variants on these cogs:

```
$METASIMPATH/scripts/CopyCogs.sh
$METASIMPATH/scripts/ReverseStrand.sh
```

Then we run some additional cleaning scripts, removing cogs that have been filtered at the median coverage step and testing that strains are indeed different:
```
$METASIMPATH/FilterTau.sh
$METASIMPATH/Degenerate.sh
```

Then we can move back to the DESMAN scg results directory and compare the predicted variants for each cluster with the references:
```
cp ClusterSpecies.txt $METASIMPATHWD/SCG_Analysis
cd $METASIMPATHWD/SCG_Analysis
sed -i '1d' ClusterSpecies.txt
$METASIMPATH/scripts/VarResults.sh > VarResults.csv
```
The above file gives the number of variants that should have been and were detected in each cluster given its mapping to a reference species:
```
NV,Cluster0,1613,0
NV,Cluster107,28025,0
Cluster108,1744,55,65,57,0.846153846154,0.964912280702
NV,Cluster109,316,3
NV,Cluster111,168695,3
Cluster112,1681,134,169,146,0.792899408284,0.917808219178
NV,Cluster116,366648,0
NV,Cluster117,1718,0
NV,Cluster12,644,0
Cluster121,2095,19,21,20,0.904761904762,0.95
```

NV in the first column indicates that no variants were expected because this was a single strain reference, then the 2nd column gives the cluster name, the third column the species id and the 4th column the number of variants actually detected.
Cluster","Strain","TP","A","P","Recall","Precn")

When strains are present the format is slightly different:
1. Cluster id
2. Species id 
3. True positives, variants correctly predicted
4. Actual no. 
4. No. predicted, TP + FP
5. Recall
6. Precision


## Get core gene sequence

Go to haplotype analysis directory:
```
cd $METASIMPATHWD
mkdir CoreFasta
cd CoreFasta
```

Then run the following bash script:
```

while IFS=, read cluster strain G H R acc prediction
do
    stub=${prediction%Filtered_Tau_star.csv}

    SEL_RUN=$METASIMPATHWD/SCG_Analysis/${cluster}_scg/${stub}
    mkdir ${cluster}
    cd ${cluster}

    cut -d"," -f 1 < $METASIMPATHWD/SCG_Analysis/${cluster}_scg/${cluster}_scgcogf.csv | sort | uniq | sed '1d' > coregenes.txt
    
    python $DESMAN/scripts/GetVariantsCore.py $METASIMPATHWD/Annotate/final_contigs_gt1000_c10K.fa $METASIMPATHWD/Split/${cluster}/${cluster}_core.cogs $SEL_RUN/Filtered_Tau_star.csv coregenes.txt -o ${
cluster} 

    cd ..

done < ../ValidateStrains.csv
```

It will generate haplotype sequences for each gene in each directory. If you need reversed cogs in the 5-3' direction for constructing phylogenies with references use the '-r' flag.




## Validate haplotype prediction

We begin by selecting clusters which are expected to have variants i.e. map to species with multiple strains:

```
grep -v "NV," VarResults.csv > VarResults_V.csv
```

```
$METASIMPATH/scripts/ValidateStrains.sh > ValidateStrains.csv
```

```
$METASIMPATH/scripts/CompareTau.sh > CompareTau.csv
```

## Assign genes to genomes

```
mkdir AssignGenes

cut -d"," -f1 < ../SCG_AnalysisF/CompareTau10.txt > ClusterV.txt

./CatGenes.sh > AllV.genes

sed 's/>\(\S\+\)\s.*$/>\1/' ../AssignContigs/AllGenomes.fasta > AllGenomesR.fasta

nohup python ./gene_read_count_per_genome.py AllV.genes AllGenomesR.fasta ../Map/*mapped.sorted.bam > AllV_counts.tsv&
```

## Detect variants on all genes


Get base frequencies at all positions, in all clusters with 75% of SCGs in single-copies:
```
mkdir AllFreq
while read -r cluster 
do
    echo $cluster

    python $DESMAN/scripts/ExtractCountFreqGenes.py -g Split/${cluster}/${cluster}.genes ./SplitBam/${cluster}/ReadcountFilter --output_file AllFreq/${cluster}.freq > AllFreq/${cluster}log.txt&

done < Concoct/Cluster75.txt
```

Run variant detection across all genes for all clusters:
```
cd AllFreq

for file in *.freq
do
    stub=${file%.freq}
    #stub=${stub#./Variants\/}
#Variants/Cluster32_scg.freq
    echo $stub
    mkdir ${stub}
    
    cp $file ${stub}
    cd ${stub}    

    python $DESMAN/desman/Variant_Filter.py ${stub}.freq -o ${stub} -m 1. -f 25.0 -p > ${stub}.vout&
    
    cd ..
done
```