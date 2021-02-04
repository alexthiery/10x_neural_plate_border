#!/bin/bash
#SBATCH --job-name=10x-NPB
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

export NXF_VER=20.07.1

nextflow run ./NF-downstream_analysis/main.nf \
--samplesheet ./samplesheets/NF-downstream_analysis/samplesheet.csv \
-profile crick \
-resume