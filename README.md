## 10x Neural Plate Border

This repository contains documentation and source code for the analysis in [Thiery et al. 2023](https://doi.org/10.7554/eLife.82717).

### scRNAseq Alignment

Alignment was performed using a custom Nextflow pipeline that tagged mitochondrial and sex chromosome genes for downstream analysis. Alignment utilized CellRanger 4.0.0, and velocity counts were generated from BAM files using Velocyto. The main workflow is available [here](NF-scRNAseq_alignment/main.nf).

### Downstream Processing Overview

The downstream analysis pipeline was built in Nextflow to allow for parallel processing and reproducibility. This pipeline involves many steps, including splitting and parallel processing of data subsets. The main Nextflow workflow can be found [here](NF-downstream_analysis/main.nf). In short, the main steps involved in the analysis are as follows:

- **Data Integration**: Combined sequencing runs using [STACAS](https://doi.org/10.1093/bioinformatics/btaa755).
  - Scripts: [Preprocessing](NF-downstream_analysis/bin/seurat/1_preprocess.R), [Integration](NF-downstream_analysis/bin/seurat/2_integration.R), [Quality Control](NF-downstream_analysis/bin/seurat/3_integration_qc.R)

- **Regression of Confounding Factors**: Controlled for sex and cell cycle effects.
  - Scripts: [Sex Filtering](NF-downstream_analysis/bin/seurat/4_sex_filt.R), [Cell Cycle Adjustment](NF-downstream_analysis/bin/seurat/5_cell_cycle.R)

- **Contamination Removal**: Excluded unwanted cell populations.
  - [Script](NF-downstream_analysis/bin/seurat/6_contamination_filt.R)

- **Data Splitting and Parallel Processing**:
  - The dataset was split by developmental timepoint, and pre-placodal and neural crest clusters were subset and reclustered for detailed neural plate border analysis. Each subset and the full dataset were processed in parallel through:
    - **Reclustering**
      - [Script](NF-downstream_analysis/bin/seurat/subset_cluster.R)
    - **Cell State Classification**
      - [Script](NF-downstream_analysis/bin/seurat/state_classification.R)
    - **[Antler](https://github.com/juliendelile/Antler) Gene Module Calculation**
      - [Script](NF-downstream_analysis/bin/other/gene_modules.R)
    - **RNA Velocity Analysis**
      - Scripts: [Seurat to h5ad Conversion](NF-downstream_analysis/modules/local/seurat_h5ad/bin/seurat_h5ad.R), [Loom Creation](NF-downstream_analysis/modules/local/scvelo/seurat_intersect_loom/bin/seurat_intersect_loom.py), [scVelo Analysis](NF-downstream_analysis/modules/local/scvelo/scvelo_run/bin/scvelo_run.py)
    - **[CellRank](https://www.nature.com/articles/s41592-021-01346-6) Cell Lineage Probability Analysis**
      - [Script](NF-downstream_analysis/modules/local/scvelo/cellrank_run/bin/cellrank_run.py)

- **Gene Module Analysis and Visualisation**:
  - Gene module heatmap and latent-time lineage dynamics for key lineage modules across the full dataset (Figure 4E).
    - [Script](NF-downstream_analysis/bin/other/gene_modules_subset_latent_time.R)
  - Gene module heatmap and latent-time lineage dynamics across the neural plate border subset (Figure 7E).
    - [Script](NF-downstream_analysis/bin/other/gene_modules_npb_latent_time.R)

- **Co-expression Analysis**:
  - Co-expression of neural plate border gene modules across various cell clusters from the full dataset (Figure 6).
    - [Script](NF-downstream_analysis/bin/other/coexpression_analysis_npb.R)

- **Spatial expression dynamics**:
  - In silico spatial expression modelling (Figure 5A/B).
    - [Script](NF-downstream_analysis/bin/other/spatial_expression_modelling.R)
  - HCR analysis (Figure 5D).
    - [Script](NF-downstream_analysis/bin/other/hcr_intensity_plot.R)