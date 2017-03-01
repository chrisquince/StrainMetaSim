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
```
