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
    ```

4.




