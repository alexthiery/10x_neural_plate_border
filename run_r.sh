#!/bin/bash
#SBATCH --job-name=otic-reprogamming
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

nextflow pull alexthiery/10x_neural_plate_border -r R

export NXF_VER=20.07.1

nextflow run alexthiery/10x_neural_plate_border \
-r R \
-profile crick