#!/bin/bash
#SBATCH --job-name=10x-NPB
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.6
ml Singularity/3.4.2
ml Graphviz

export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/singularity

nextflow run ./NF-downstream_analysis/main.nf \
--input ./NF-downstream_analysis/samplesheet.csv \
--outdir ./output/NF-downstream_analysis \
--debug \
-profile crick \
-resume