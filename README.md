# StrainMetaSim

This repository contains scripts and downloads to data necessary to recreate the `complex strain' in silico mock community from Quince et al "DESMAN: a new tool to decipher strain-level divergence de novo from metagenomes".

1. Lets create a directory for working in and set PATH to scripts:
    ```
        mkdir ComplexStrainSim
        cd ComplexStrainSim
        export METASIMPATH=/MyPath/StrainMetaSim/scripts
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




