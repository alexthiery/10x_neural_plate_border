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

nextflow run ./NF-scRNAseq_alignment/main.nf \
--input ./NF-scRNAseq_alignment/crick_samplesheet.csv \
--outdir ./output/NF-scRNAseq_alignment \
-profile crick \
-resume