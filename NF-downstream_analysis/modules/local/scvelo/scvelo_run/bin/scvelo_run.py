#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import scvelo as scv
import cellrank as cr
import numpy as np
import warnings
import scanpy as sc
import pandas as pd

def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help="Input loom path.", metavar='')
    parser.add_argument('-o', '--output', type=str, help="Prefix for output files.", metavar='')
    parser.add_argument('-m', '--velocityMode', type=str, help="Method for calculating RNA velocity. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.", default='dynamical', metavar='')
    parser.add_argument('-c', '--clusterColumn', type=str, help="Name of cluster column.", default='seurat_clusters', metavar='')
    parser.add_argument('-s', '--stageColumn', type=str, help="Name of stage column.", default=None, metavar='')
    parser.add_argument('-b', '--batchColumn', type=str, help="Name of batch column.", default=None, metavar='')
    parser.add_argument('-g', '--genes', help="Genes of interest to plot on velocity.", nargs='+')
    parser.add_argument('-n', '--ncores', type=int, help="Number of cores used for parallelisation.", metavar='', default=1)
    parser.add_argument('-d', '--dpi', type=int, help='Set DPI for plots.', default='240')
    parser.add_argument('-gd', '--geneDiscovery', type=bool, help='Run unbiased gene discovery.', default=True)
    parser.add_argument('-k', '--diffKinetics', type=bool, help='Whether to run diff kinetics.', default=False)
    parser.add_argument('-f', '--forceCol', type=bool, help='Force plotting color bars.', default=False)
    parser.add_argument('-cc', '--coloursColumn', type=str, help='Name of cell colours column.', default=None)
    parser.add_argument('-np', '--npcs', type=int, help='Number of PCs to use for calculating moments', default=30)
    parser.add_argument('-nn', '--nneighbours', type=int, help='Number of neighbours to use for calculating moments', default=30)
    parser.add_argument('-r', '--root', type=str, help='Name of root', default=None)
    parser.add_argument('-re', '--rootEarliest', type=str, help='Space separated array specifying temporal arrangement of stages (i.e. hh4,hh5,hh6)', nargs='+', default=None)
    parser.add_argument('-rc', '--rootCol', type=str, help='Name of root metadata column', default=None)
    parser.add_argument('-e', '--end', type=str, help='Name of end', default=None)
    parser.add_argument('-el', '--endLatest', type=str, help='Space separated array specifying temporal arrangement of stages (i.e. hh4,hh5,hh6)', nargs='+', default=None)
    parser.add_argument('-ec', '--endCol', type=str, help='Name of end metadata column', default=None)
    parser.add_argument('-w', '--weightDiffusion', type=float, help='Weight applied to couple latent time with diffusion-based velocity pseudotime', default=None)
    return parser.parse_args(args)

def check_args(args, adata):
    if not os.path.isfile(args.input):
        raise Exception(f"'--input': '{args.input}' is not a valid path")
        
    if args.velocityMode not in ['dynamical', 'deterministic', 'stochastic']:
        raise Exception(f"'--velocityMode' must be set to either: 'dynamical', 'deterministic', or 'stochastic'")
    
    for arg in [args.clusterColumn, args.stageColumn, args.batchColumn, args.coloursColumn, args.rootCol]:
        if arg is not None:
            if arg not in adata.obs.columns:
                raise Exception(f"'{args.input}' is not a column in adata.obs")
                
    if args.root is not None and args.rootCol is None:
        raise Exception(f"'--rootCol' must be set when '--root' is specified")
    
    if args.rootEarliest is not None:
        if args.rootCol is None:
            raise Exception(f"'--rootCol' must be set when '--rootEarliest' is specified")
        if args.root is not None:
            warnings.warn("both '--root' and '--rootEarliest' have been set. '--rootEarliest' is prioritised.")
        args.root = [stage for stage in args.rootEarliest if stage in adata.obs.stage.unique()][0]

    if args.end is not None and args.endCol is None:
        raise Exception(f"'--endCol' must be set when '--end' is specified")
        
    if args.endLatest is not None:
        if args.endCol is None:
            raise Exception(f"'--endCol' must be set when '--endLatest' is specified")
        if args.end is not None:
            warnings.warn("both '--end' and '--endLatest' have been set. '--endLatest' is prioritised.")
        args.end = [stage for stage in args.endLatest if stage in adata.obs.stage.unique()][len(adata.obs.stage.unique()) -1]
        
# Read in loom data
def read_loom(loom_path, clusterColumn, stageColumn, batchColumn):
    adata = scv.read(loom_path)
    # set var columns as categorical variables
    adata.obs[clusterColumn]=adata.obs[clusterColumn].astype('category')
    if stageColumn is not None:
        adata.obs[stageColumn]=adata.obs[stageColumn].astype('category')
    if batchColumn is not None:   
        adata.obs[batchColumn]=adata.obs[batchColumn].astype('category')
    return(adata)

# extract unique elements in a list whilst preserving their order of appearance in the list
def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

# Plot splice proportions
def plot_proportions(adata, clusterColumn):
    ax = scv.pl.proportions(adata, groupby=clusterColumn, show=False, save='proportions.pdf')

# scVelo pre-processing
def preprocess_anndata(adata):
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)
    return(adata)

# calculate cell velocity
def calc_velocity(adata, velocityMode, ncores, groupby=None, diffKinetics=False):
    if velocityMode not in ['dynamical', 'deterministic', 'stochastic']:
        Exception(f"'--velocityMode': '{velocityMode}' is not valid. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.")
    if velocityMode == 'dynamical' and diffKinetics == False:
        scv.tl.recover_dynamics(adata, n_jobs=ncores)
    scv.tl.velocity(adata, mode=velocityMode, groupby=groupby, diff_kinetics=diffKinetics)
    scv.tl.velocity_graph(adata, n_jobs=ncores)
    return(adata)

def plot_velocity(adata, clusterColumn, threshold=.1, arrow_length=5, arrow_size=2, plot_dir="", prefix="", dpi=240):
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    scv.pl.velocity_embedding(adata, color=clusterColumn, arrow_length=arrow_length, arrow_size=arrow_size, save=prefix+'velocity_embedding.png', dpi=dpi)
    scv.pl.velocity_embedding_grid(adata, color=clusterColumn, arrow_length=arrow_length, arrow_size=arrow_size, basis='umap', save=prefix+'velocity_embedding_grid.png', dpi=dpi)
    scv.pl.velocity_embedding_stream(adata, color=clusterColumn, basis='umap', save=prefix+'velocity_embedding_stream.png', dpi=dpi)
    scv.pl.velocity_graph(adata, threshold=threshold, color=clusterColumn, basis='umap', save=prefix+'velocity_graph.png', dpi=dpi)
    scv.settings.figdir = keep_figdir
    
# plot expression of genes of interest across velocity UMAPs
def plot_genes(adata, genes_dict, clusterColumn, dpi=240, plot_dir="", prefix=""): # genes_dict is a key value pair with {name:gene(s)}
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    for key, value in genes_dict.items():
        scv.pl.velocity(adata, value, save=prefix+key+'.png', color=clusterColumn, dpi=dpi)
    scv.settings.figdir = keep_figdir

# plot expression of genes of interest across velocity UMAPs
def plot_genes_dynamical(adata, genes_dict, clusterColumn, dpi=240, plot_dir="", prefix=""): # genes_dict is a key value pair with {name:gene(s)}
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    
    if plot_dir != "" and not os.path.exists(scv.settings.figdir+plot_dir):
        os.makedirs(scv.settings.figdir+plot_dir)
    for key, value in genes_dict.items():
        scv.pl.scatter(adata, value, ylabel=key, color=clusterColumn, frameon=False, save=prefix+key+'_dynamical.png', dpi=dpi)
        scv.pl.scatter(adata, x='latent_time', y=value, color=clusterColumn, save=prefix+key+'_dynamical_latent_time.png', frameon=False)

        
# identify genes which explain the vector field in a given lineage
def ident_top_genes(adata, clusterColumn, min_corr=.3, n_genes=5):
    scv.tl.rank_velocity_genes(adata, groupby=clusterColumn, min_corr=min_corr)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    
    top_cluster_genes = {'Cluster-'+column: df.loc[:5,column].values for column in df}
    return(top_cluster_genes) # return dictionary: {cluster_1:top_genes, cluster_2:top_genes}

def ident_top_genes_dynamical(adata, clusterColumn, n_genes=5):
    scv.tl.rank_dynamical_genes(adata, groupby=clusterColumn)
    df = scv.get_df(adata, 'rank_dynamical_genes/names')
    
    top_cluster_genes = {'Cluster-'+column: df.loc[:5,column].values for column in df}
    return(top_cluster_genes) # return dictionary: {cluster_1:top_genes, cluster_2:top_genes}

def plot_differentiation(adata, clusterColumn, dpi=240):
    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save='differentiation_scatter.png', dpi=dpi)

    
# def velocity_pseudotime(adata, dpi=240):
#     scv.tl.velocity_pseudotime(adata)
#     scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='pseudotime.png', dpi=dpi)
#     return(adata)


def paga(adata, clusterColumn, time_prior, dpi=240):
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    scv.tl.paga(adata, groups=clusterColumn)
    paga_df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    paga_df.style.background_gradient(cmap='Blues').format('{:.2g}')
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, use_time_prior='latent_time', save='paga.png', dpi=dpi)
    return(adata, paga_df)

def latent_time(adata, args):
    if args.end is not None:
        adata.obs['end_val'] = pd.Categorical(np.where(adata.obs[args.endCol] == args.end, '0', None))
        end_key = 'end_val'
    else:
        end_key = None
        
    if args.root is not None:
        adata.obs['root_val'] = pd.Categorical(np.where(adata.obs[args.rootCol] == args.root, '0', None))
        scv.tl.latent_time(adata, root_key='root_val', end_key=end_key, weight_diffusion = args.weightDiffusion)
    else:
        scv.tl.latent_time(adata, end_key=end_key, weight_diffusion = args.weightDiffusion)
        
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png', dpi=args.dpi)
    return(adata)


def plot_latent_time_heatmap(adata, clusterColumn, stageColumn=None, batchColumn=None, forceCol=False, plot_dir=""):
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    col_color=[clusterColumn, stageColumn, batchColumn]
    # set annotation columns based on presence of multiple levels in cat var
    if forceCol == True:
        col_color = [i for i in col_color if i]
    else:
        col_color = [i for i in col_color if i if len(set(adata.obs[i])) >1]

    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color=col_color, n_convolve=100, save='latent_time.png')
    scv.settings.figdir = keep_figdir


def calc_diff_kinetics(adata, clusterColumn, plot_dir=""):
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
    scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby=clusterColumn)
    scv.pl.scatter(adata, basis=top_genes[:20], ncols=5, add_outline='fit_diff_kinetics', frameon=False, add_linfit=True, linewidth=2, save="top_genes_diff_kinetics_1-20.png")
    scv.pl.scatter(adata, basis=top_genes[20:40], ncols=5, add_outline='fit_diff_kinetics', frameon=False, add_linfit=True, linewidth=2, save="top_genes_diff_kinetics_21-40.png")
    scv.settings.figdir = keep_figdir
    return(adata)


def directed_paga(adata, clusterColumn, plot_dir="", prefix="", dpi=240, weight_connectivities=0.2):
    keep_figdir = scv.settings.figdir # Save original plot path
    scv.settings.figdir = scv.settings.figdir+plot_dir
    scv.tl.paga(adata, groups=clusterColumn, root_key="initial_states_probs", end_key="terminal_states_probs",
                use_time_prior="velocity_pseudotime")
    cr.pl.cluster_fates(adata, mode="paga_pie", cluster_key=clusterColumn, basis="umap",
                        legend_kwargs={"loc": "top right out"}, legend_loc="top left out", node_size_scale=5,
                        edge_width_scale=1, max_edge_width=4, title="directed PAGA",
                        save=prefix+'directed_paga.png', dpi=dpi)
    scv.settings.figdir = keep_figdir

    
# functions for running scvelo/cellrank subworkflows
def run_scvelo_deterministic(adata, args):
    if args.geneDiscovery == True:
        # Identify and plot genes that have cluster-specific differential velocity expression
        top_cluster_genes = ident_top_genes(adata, args.clusterColumn)
        plot_genes(adata=adata, genes_dict=top_cluster_genes, clusterColumn=args.clusterColumn, plot_dir = 'gene_discovery/', dpi=args.dpi)
        plot_differentiation(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
        
    # Plot velocity for manual GOI
    if args.genes is not None:
        manual_genes = {args.genes[i]: args.genes[i] for i in range(0, len(args.genes))}
        plot_genes(adata=adata, genes_dict=manual_genes, clusterColumn=args.clusterColumn, plot_dir = 'goi/', dpi=args.dpi)
    
    return(adata)


def run_scvelo_dynamical(adata, args):
    adata = latent_time(adata, args)
    if args.geneDiscovery == True:
        plot_latent_time_heatmap(adata, clusterColumn = args.clusterColumn, stageColumn = args.stageColumn, forceCol = args.forceCol, batchColumn = args.batchColumn, plot_dir='gene_discovery_dynamical/')
        # Identify and plot genes that have cluster-specific differential velocity expression
        top_cluster_genes = ident_top_genes_dynamical(adata, clusterColumn = args.clusterColumn)
        plot_genes_dynamical(adata, genes_dict=top_cluster_genes, clusterColumn=args.clusterColumn, plot_dir = 'gene_discovery_dynamical/', dpi=args.dpi)
    # Plot velocity for manual GOI
    if args.genes is not None:
        manual_genes = {args.genes[i]: args.genes[i] for i in range(0, len(args.genes))}
        plot_genes_dynamical(adata=adata, genes_dict=manual_genes, clusterColumn=args.clusterColumn, plot_dir = 'goi_dynamical/', dpi=args.dpi)

    print('Calculating differential kinetics and recomputing velocity')
    adata = calc_diff_kinetics(adata=adata, clusterColumn=args.clusterColumn, plot_dir="")
    adata = calc_velocity(adata=adata, velocityMode=args.velocityMode, ncores=args.ncores, diffKinetics=True, groupby=args.clusterColumn)
    plot_velocity(adata=adata, clusterColumn=args.clusterColumn, plot_dir="diff_kinetics/", dpi=args.dpi)
    return(adata)

def main(args=None):
    
    args = parse_args(args)
    
    # Set global settings
    scv.logging.print_version()
    scv.settings.plot_prefix = "" # remove plot prefix
    scv.settings.n_jobs = args.ncores  # set max width size for presenter view
    
    adata = read_loom(loom_path=args.input, clusterColumn=args.clusterColumn, stageColumn=args.stageColumn, batchColumn=args.batchColumn)
    
    check_args(args, adata)
    
    # Set plotting colours if available
    if args.coloursColumn is not None:
        print('Setting cluster colours using ' + args.coloursColumn + ' column')
        # Extract cell colours for plotting
        adata.uns[args.clusterColumn + '_colors'] = unique(adata.obs.sort_values(by=[args.clusterColumn], inplace=False)[args.coloursColumn])
    
    plot_proportions(adata=adata, clusterColumn=args.clusterColumn)
    
    adata = preprocess_anndata(adata=adata)
    # calculate means and uncentered variances across NN in PCA space
    scv.pp.moments(data=adata,  n_pcs=args.npcs, n_neighbors=args.nneighbours)
    adata = calc_velocity(adata=adata, velocityMode=args.velocityMode, ncores=args.ncores)
    plot_velocity(adata=adata, clusterColumn=args.clusterColumn, dpi=args.dpi)

    # Run dynamical or deterministic models and plot genes
    if args.velocityMode in ['deterministic', 'stochastic']:
        adata = run_scvelo_deterministic(adata, args)
        adata, paga_df = paga(adata, clusterColumn=args.clusterColumn, time_prior='velocity_pseudotime', dpi=args.dpi)
    
    if args.velocityMode == 'dynamical':
        adata = run_scvelo_dynamical(adata, args)
        adata, paga_df = paga(adata, clusterColumn=args.clusterColumn, time_prior='latent_time', dpi=args.dpi)

    adata.write(args.output + '_scvelo.h5ad')
    # return(args, adata)

if __name__ == '__main__':
    sys.exit(main())



# args = ['-i', '../output/NF-downstream_analysis_stacas/scvelo/NF-scRNAseq_alignment_out/seurat_intersect_loom/NF-scRNAseq_alignment_out_seurat_intersect.loom', '-o', 'out.h5ad', '-m', 'dynamical', '-c',
#         'scHelper_cell_type', '-s', 'stage', '-b', 'run', '--ncores', '32', '--coloursColumn', 'cell_colours', '--npcs', '20', '--nneighbours', '20', '--root', 'hh4', '--rootCol', 'stage',
#         '--weightDiffusion', '0.2', '--diffKinetics', 'True']

# args, adata = main(args)