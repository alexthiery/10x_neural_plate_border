#!/usr/bin/env python

import os
import sys
import errno
import argparse
import loompy
import scvelo as scv
import pandas as pd


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    parser.add_argument('-si', '--seuratIntersect', help="Pre-filtered H5 seurat object for intersecting with loom", metavar='')
    parser.add_argument('-a', '--annotations', help="Pre-filtered H5 seurat object for intersecting with loom", metavar='')
    parser.add_argument('-lo', '--loomOutput', help="Path to concatenated loom output file.", metavar='')
    return parser.parse_args(args)
    

def check_args(args=None):
    if os.path.isdir(args.input):
        if args.loomOutput is None:
            raise Exception("'--loomOutput' must be specified when '--input' is a directory containing loom files to concatenate.")
    elif not os.path.isfile(args.input):
        raise Exception(f"'--input': '{args.input}' is not a valid path.")


def get_file_paths(path):
    files = [path+'/'+f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    return files


# Subset loom files from directory and concatenate if necessary
def concatenate_loom(args=None):
    input_path = os.path.abspath(args.input)

    if os.path.isdir(input_path):
        input_files = get_file_paths(input_path)
        selected_files = [file for file in input_files if file.lower().endswith('loom')]

        if len(selected_files) < len(input_files):
            if(len(selected_files) == 0):
                raise Exception("'--input' must contain '.loom' files.")
            print("\n'--input' directory contains non .loom files.\nThe following files have been subset for merging:", selected_files)

        output_file = os.path.abspath(args.loomOutput)

        print("\nCombining loom files:", selected_files, "\nOutput loom file:", output_file, "\n")
        loompy.combine(selected_files, output_file=output_file, key="Accession")
        return(output_file)

    elif os.path.isfile(input_path):
        return(input_path)


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


def preprocess_anndata(adata):
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

def calc_moments(adata):
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000, enforce=True)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


def main(args=None):
    args = parse_args(args)    
    check_args(args)
    
    # Generate merged loom file if required
    loom_path = concatenate_loom(args)

    # Read in loom data
    adata = read_loom(loom_path)

    # If seurat and annotations are provided, then subset adata names and cells accordingly and merge
    if args.seuratIntersect:
        merge_data = merge_seurat_loom(adata, args.annotations, args.seuratIntersect)

    # scVelo pre-processing
    preprocess_anndata(merge_data)

    calc_moments(merge_data)

    print(merge_data)

if __name__ == "__main__":
    sys.exit(main())
