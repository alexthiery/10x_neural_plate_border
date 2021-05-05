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
    parser.add_argument('-vm', '--velocityMode', help="Method for calculating RNA velocity. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.", default='dynamical', metavar='')
    parser.add_argument('-o', '--output', help="Path to concatenated loom output file.", metavar='')
    return parser.parse_args(args)
    

# def check_args(args=None):
#     if os.path.isdir(args.input):
#         if args.loomOutput is None:
#             raise Exception("'--loomOutput' must be specified when '--input' is a directory containing loom files to concatenate.")
#     elif not os.path.isfile(args.input):
#         raise Exception(f"'--input': '{args.input}' is not a valid path.")

def read_loom(loom_path):
    adata = scv.read(loom_path)
    return(adata)

def preprocess_anndata(adata):
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

def calc_moments(adata):
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

def calc_velocity(args, adata):
    if args.velocityMode == 'dynamical':
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
    elif args.velocityMode == 'deterministic':
        scv.tl.velocity(adata, mode='deterministic')
    elif args.velocityMode == 'stochastic':
        scv.tl.velocity(adata)
    else:
        Exception(f"'--velocityMode': '{args.velocityMode}' is not valid. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.")

    scv.tl.velocity_graph(adata)

def plot_velocity(adata):
    scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, save=True)
    scv.pl.velocity_embedding_grid(adata, basis='umap', save=True)
    scv.pl.velocity_embedding_stream(adata, basis='umap', save=True)

def latent_time(adata):
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save=True)
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100)


def main(args=None):
    args = parse_args(args)
    # check_args(args)

    # Read in loom data
    adata = read_loom(args.input)

    # scVelo pre-processing
    preprocess_anndata(adata)

    calc_moments(adata)


    calc_velocity(args, adata)

    plot_velocity(adata)

    latent_time(adata)

if __name__ == "__main__":
    sys.exit(main())
