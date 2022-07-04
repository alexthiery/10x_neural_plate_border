#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import warnings

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)

def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, help="Input hdf5 path.", metavar='')
    parser.add_argument('-o', '--output', type=str, help="Prefix for output files.", metavar='')
    parser.add_argument('-c', '--clusterColumn', type=str, help="Name of cluster column.", default='seurat_clusters', metavar='')
    parser.add_argument('-n', '--ncores', type=int, help="Number of cores used for parallelisation.", metavar='', default=1)
    parser.add_argument('-d', '--dpi', type=int, help='Set DPI for plots.', default='240')
    parser.add_argument('-ck', '--combineKernel', type=bool, help='Combine velocity and connectivity transition matrices', default=True)
    parser.add_argument('-kp', '--kernelProportion', type=float, help='The proportional weighting of the connectivity kernel in the combined transition matrix', default=0.2)
    parser.add_argument('-sm', '--scaleMax', type=float, help='Softscale max value for when calculating transition matrix', default=None)
    parser.add_argument('-dt', '--dataType', type=str, help='Which data subset the analysis is carried out on', default='full')
    return parser.parse_args(args)

# Read in h5ad data
def read_h5ad(h5ad_path):
#     adata = scv.read(h5ad_path, X_name='')
    adata = scv.read(h5ad_path)
    return(adata)
    
def write_lineage_probs(adata):
    for lineage in adata.obsm['to_terminal_states'].names:
        colname = 'lineage_' + lineage + '_probability'
        if colname in adata.obs.columns:
            warnings.warn(colname + ' is already specified in adata.obs. Overriding original entry.')
            del adata.obs[colname]
        
        adata.obs[colname] = adata.obsm['to_terminal_states'][lineage]
    
    print('Adding terminal states to metadata')
    adata.obs = adata.obs.reindex(copy=False)
    return(adata)

def allDataTerminalStates(adata, estimator, dpi):
    estimator.set_terminal_states({"neural": adata[adata.obs["scHelper_cell_type"].isin(['HB', 'MB', "FB"]) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names,
                  "NC": adata[adata.obs["scHelper_cell_type"].isin(['dNC', 'NC']) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names,
                  "placodal": adata[adata.obs["scHelper_cell_type"].isin(['aPPR', 'pPPR', 'PPR']) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names})
    cr.pl.terminal_states(adata, save='terminal_states.pdf', dpi=dpi)
    return(estimator)

def transferLabelTerminalStates(adata, estimator, dpi):
    estimator.set_terminal_states({"neural": adata[adata.obs["scHelper_cell_type"].isin(['HB', 'MB', "FB"]) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names,
                  "NC": adata[adata.obs["scHelper_cell_type"].isin(['dNC', 'NC']) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names,
                  "placodal": adata[adata.obs["scHelper_cell_type"].isin(['aPPR', 'pPPR', 'PPR']) & adata.obs["stage"].isin(['ss8', 'ss4'])].obs_names})
    cr.pl.terminal_states(adata, save='terminal_states.pdf', dpi=dpi)
    return(estimator)

def terminalStates_ppr_nc(adata, estimator, dpi):
    estimator.set_terminal_states({"NC": adata[adata.obs["scHelper_cell_type"].isin(['dNC', 'NC']) & adata.obs["stage"].isin(['ss8'])].obs_names,
                  "placodal": adata[adata.obs["scHelper_cell_type"].isin(['aPPR', 'pPPR', 'PPR']) & adata.obs["stage"].isin(['ss8'])].obs_names})
    cr.pl.terminal_states(adata, save='terminal_states.pdf', dpi=dpi)
    return(estimator)

def write_lineage_probs(adata):
    for lineage in adata.obsm['to_terminal_states'].names:
        colname = 'lineage_' + lineage + '_probability'
        if colname in adata.obs.columns:
            warnings.warn(colname + ' is already specified in adata.obs. Overriding original entry.')
            del adata.obs[colname]
        
        adata.obs[colname] = adata.obsm['to_terminal_states'][lineage]
    
    print('Adding terminal states to metadata')
    adata.obs = adata.obs.reindex(copy=False)
    return(adata)


def main(args=None):
    args = parse_args(args)
    # check_args(args)


    # Set global settings
    scv.logging.print_version()
    scv.settings.plot_prefix = "" # remove plot prefix
    scv.settings.n_jobs = args.ncores  # set max width size for presenter view

    adata = read_h5ad(h5ad_path=args.input)

    if args.combineKernel == True:
        combined_kernel = cr.tl.transition_matrix(adata, weight_connectivities=args.kernelProportion, softmax_scale=args.scaleMax, show_progress_bar=False)
        g = cr.tl.estimators.GPCCA(combined_kernel)

    if args.dataType == 'full':
        g = allDataTerminalStates(adata, g, dpi=args.dpi)
    elif args.dataType == 'labelTransfer':
        g = transferLabelTerminalStates(adata, g, dpi = args.dpi)
    elif args.dataType == 'ppr_nc':
        g = terminalStates_ppr_nc(adata, g, dpi=args.dpi)
        
    g.compute_absorption_probabilities()
    cr.pl.lineages(adata, same_plot=False, save='absorption_probabilities.pdf', dpi = args.dpi)
    
    if all(value is None for value in adata.obs['terminal_states_probs']):
        adata.obs = adata.obs.drop(['terminal_states_probs'], axis=1)

    adata = write_lineage_probs(adata)
    adata.write(args.output + '_cellrank.h5ad')
    adata.obs.to_csv(args.output + '_metadata.csv')
    # return(args, adata)

if __name__ == '__main__':
    sys.exit(main())

# # set args for interactive testing
# args = ['-i', '../output/NF-downstream_analysis_stacas/scvelo/NF-scRNAseq_alignment_out/scvelo_run/NF-scRNAseq_alignment_out_scvelo.h5ad', '-o', 'out.h5ad', '-c',
#         'seurat_clusters', '--ncores', '8', '-ck', 'True', '-kp', '0.8']