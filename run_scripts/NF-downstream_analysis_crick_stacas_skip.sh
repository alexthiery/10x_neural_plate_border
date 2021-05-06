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
--input ./NF-downstream_analysis/samplesheet.csv \
--outdir ./output/NF-downstream_analysis_stacas \
--debug \
--integration STACAS \
--skip_seurat_filtering \
--seurat_out /camp/home/thierya/scratch/10x_neural_plate_border/output/NF-downstream_analysis_stacas/seurat/6_contamination_filt/rds_files/contamination_filt_data.RDS \
--seurat_annotations /camp/home/thierya/scratch/10x_neural_plate_border/output/NF-downstream_analysis_stacas/seurat/1_integration/seurat_annotations.csv \
-profile crick \
-resume