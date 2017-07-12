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

## Running DESMAN on the complex mock

This assumes that DESMAN and CONCOCT are installed and their paths 
set to the variables DESMAN and CONCOCT respectively e.g. (changing paths to your 
system set-up):

```
export CONCOCT=~/Installed/CONCOCT/
export DESMAN=/mnt/gpfs/chris/repos/DESMAN/
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
concoct_refine clustering_gt1000.csv original_data_gt1000.csv clustering_gt1000_scg_sort.csv > concoct_ref.out
```

This should result in XXX clusters with 75% single copy copy SCGs:
```
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c cluster_refine.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv
```

