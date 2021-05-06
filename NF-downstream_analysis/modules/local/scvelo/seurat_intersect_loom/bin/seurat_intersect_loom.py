#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import loompy
import scvelo as scv
import pandas as pd


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-li', '--loomInput', help="Input loom file.", metavar='')
    parser.add_argument('-si', '--seuratInput', help="Input loom file.", metavar='')
    parser.add_argument('-a', '--annotations', help="Input loom file.", metavar='')
    parser.add_argument('-o', '--output', help="Input loom file.", metavar='')
    return parser.parse_args(args)
    

# def check_args(args=None):
#     if not os.path.isdir(args.input):
#         raise Exception("'--input' path is not a directory.")

def read_loom(loom_path):
    adata = scv.read(loom_path)
    adata.var = adata.var.set_index('Accession')
    return(adata)

# subset annotations data frame by remaining genes in seurat object
def seurat_filter_annotations(seurat_hd5, annotations):
    seurat_hd5 = scv.read(seurat_hd5)
    annotations = pd.read_csv(annotations)
    annotations = annotations.loc[annotations.Gene.isin(seurat_hd5.var.features)]
    annotations = annotations.set_index('Accession')
    return(annotations, seurat_hd5)

# subset adata by annotations and set gene names
def annotation_adata_subset(adata, annotations):
    if not all(annotations.index.isin(adata.var.index)):
        raise Exception('Seurat genes missing from loom file. Check Velocyto output.')

    adata = adata[:, annotations.index]
    adata.var = pd.merge(adata.var, annotations, on='Accession')
    adata.var = adata.var.set_index('Gene')
    return(adata)

def rename_loom_cells(adata, seurat):
    adata.obs.index = adata.obs.index.str.replace("cellranger:", "")
    adata.obs.index = adata.obs.index.str[:-1]
    seurat.obs.index = seurat.obs.index.str[:-2]

    # check all seurat cell names are in velocyto object and subset adata cells by cells remaining in seurat object
    if not all(adata.obs.index.isin(seurat.obs.index)):
        # raise Exception("Seurat cell names are missing from loom object.")
        print("\nSeurat cell names are missing from loom object - please check that correct loom files are being processed.\nMatching cells have been subset for downstream analysis.")

    return(adata, seurat)


def merge_seurat_loom(loom, annotations, seurat_hd5):
    annotations, seurat = seurat_filter_annotations(seurat_hd5, annotations)
    adata = annotation_adata_subset(loom, annotations)
    adata, seurat = rename_loom_cells(adata, seurat)
    merge_data = scv.utils.merge(adata, seurat)
    return(merge_data)


def main(args=None):
    args = parse_args(args)
    # check_args(args)
    
    # Read in loom data
    adata = read_loom(args.loomInput)

    merge_data = merge_seurat_loom(adata, args.annotations, args.seuratInput)

    merge_data.write_loom(args.output, write_obsm_varm=True)

if __name__ == "__main__":
    sys.exit(main())
