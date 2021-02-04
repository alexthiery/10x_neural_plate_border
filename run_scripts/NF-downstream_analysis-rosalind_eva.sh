#!/bin/bash
#SBATCH --job-name=otic-reprogamming
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk


## LOAD REQUIRED MODULES
module load apps/singularity/3.5.3


export NXF_VER=20.07.1

nextflow run ./NF-downstream_analysis/main.nf \
    -profile rosalind_eva \
    --input ./NF-downstream_analysis/local_samplesheet.csv