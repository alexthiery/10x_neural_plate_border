#!/bin/bash
#SBATCH --job-name=10x-NPB
#SBATCH -t 72:00:00
#SBATCH -p brc
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk


## LOAD REQUIRED MODULES
module load apps/singularity/3.5.3


export NXF_VER=20.07.1

nextflow run ./NF-downstream_analysis/main.nf \
    -profile rosalind_alex \
    --input ./NF-downstream_analysis/local_samplesheet.csv