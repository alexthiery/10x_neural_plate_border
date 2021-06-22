#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import scvelo as scv


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
    # parser.add_argument('-r', '--ram', type=int, help="Max memory usage in Gigabyte.", metavar='', default=12)
    parser.add_argument('-d', '--dpi', type=int, help='Set DPI for plots.', default='240')
    parser.add_argument('-gd', '--geneDiscovery', type=bool, help='Run unbiased gene discovery.', default=True)
    parser.add_argument('-f', '--forceCol', type=bool, help='Force plotting color bars.', default=False)
    return parser.parse_args(args)

# def check_args(args=None):
#     if os.path.isdir(args.input):
#         if args.loomOutput is None:
#             raise Exception("'--loomOutput' must be specified when '--input' is a directory containing loom files to concatenate.")
#     elif not os.path.isfile(args.input):
#         raise Exception(f"'--input': '{args.input}' is not a valid path.")


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

# Plot splice proportions
def plot_proportions(adata, clusterColumn):
    ax = scv.pl.proportions(adata, groupby=clusterColumn, show=False)
    scv.pl.utils.savefig_or_show(save='proportions.pdf')

# scVelo pre-processing
def preprocess_anndata(adata):
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)
    return(adata)

# calculate means and uncentered variances across NN in PCA space
def calc_moments(adata):
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    return(adata)

# calculate cell velocity
def calc_velocity(adata, velocityMode, ncores, groupby=None, diffKinetics=False):
    if velocityMode not in ['dynamical', 'deterministic', 'stochastic']:
        Exception(f"'--velocityMode': '{velocityMode}' is not valid. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.")
    if velocityMode == 'dynamical' and diffKinetics == False:
        scv.tl.recover_dynamics(adata, n_jobs=ncores)
    scv.tl.velocity(adata, mode=velocityMode, groupby=groupby, diff_kinetics=diffKinetics)
    scv.tl.velocity_graph(adata)
    return(adata)

def plot_velocity(adata, clusterColumn, threshold=.1, arrow_length=5, arrow_size=2, plot_dir="", prefix="", dpi=240):
    if plot_dir != "" and not os.path.exists(scv.settings.figdir+plot_dir):
        os.makedirs(scv.settings.figdir+plot_dir)
    scv.pl.velocity_embedding(adata, color=clusterColumn, arrow_length=arrow_length, arrow_size=arrow_size, save=plot_dir+prefix+'velocity_embedding.png', dpi=dpi)
    scv.pl.velocity_embedding_grid(adata, color=clusterColumn, arrow_length=arrow_length, arrow_size=arrow_size, basis='umap', save=plot_dir+prefix+'velocity_embedding_grid.png', dpi=dpi)
    scv.pl.velocity_embedding_stream(adata, color=clusterColumn, basis='umap', save=plot_dir+prefix+'velocity_embedding_stream.png', dpi=dpi)
    scv.pl.velocity_graph(adata, threshold=threshold, color=clusterColumn, basis='umap', save=plot_dir+prefix+'velocity_graph.png', dpi=dpi)

# plot expression of genes of interest across velocity UMAPs
def plot_genes(adata, genes_dict, clusterColumn, dpi=240, plot_dir="", prefix=""): # genes_dict is a key value pair with {name:gene(s)}
    if plot_dir != "" and not os.path.exists(scv.settings.figdir+plot_dir):
        os.makedirs(scv.settings.figdir+plot_dir)
    for key, value in genes_dict.items():
        scv.pl.velocity(adata, value, save=plot_dir+prefix+key+'.png', color=clusterColumn, dpi=dpi)

# plot expression of genes of interest across velocity UMAPs
def plot_genes_dynamical(adata, genes_dict, clusterColumn, dpi=240, plot_dir="", prefix=""): # genes_dict is a key value pair with {name:gene(s)}
    if plot_dir != "" and not os.path.exists(scv.settings.figdir+plot_dir):
        os.makedirs(scv.settings.figdir+plot_dir)
    for key, value in genes_dict.items():
        scv.pl.scatter(adata, value, ylabel=key, color=clusterColumn, frameon=False, save=plot_dir+prefix+key+'_dynamical.png', dpi=dpi)
        scv.pl.scatter(adata, x='latent_time', y=value, color=clusterColumn, save=plot_dir+prefix+key+'_dynamical_latent_time.png', frameon=False)

        
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

    
def velocity_pseudotime(adata, dpi=240):
    scv.tl.velocity_pseudotime(adata)
    scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='pseudotime.png', dpi=dpi)
    return(adata)


def paga(adata, clusterColumn, dpi=240):
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    scv.tl.paga(adata, groups=clusterColumn)
    paga_df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    paga_df.style.background_gradient(cmap='Blues').format('{:.2g}')
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save='paga.png', dpi=dpi)
    return(adata, paga_df)

def latent_time(adata, dpi=240):
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png', dpi=dpi)
    return(adata)


def plot_latent_time_heatmap(adata, clusterColumn, stageColumn=None, batchColumn=None, forceCol=False, plot_dir=""):
    keep_figdir = scv.settings.figdir
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
    if plot_dir != "" and not os.path.exists(scv.settings.figdir+plot_dir):
        os.makedirs(scv.settings.figdir+plot_dir)
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
    scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby=clusterColumn)
    scv.pl.scatter(adata, basis=top_genes[:20], ncols=5, add_outline='fit_diff_kinetics', frameon=False, add_linfit=True, linewidth=2, save=plot_dir+"/top_genes_diff_kinetics_1-20.png")
    scv.pl.scatter(adata, basis=top_genes[20:40], ncols=5, add_outline='fit_diff_kinetics', frameon=False, add_linfit=True, linewidth=2, save=plot_dir+"/top_genes_diff_kinetics_21-40.png")
    return(adata)


def main(args=None):
    
    args = parse_args(args)
    # check_args(args)
    
    # Set global settings
    scv.logging.print_version()
    scv.settings.plot_prefix = "" # remove plot prefix
    scv.settings.n_jobs = args.ncores  # set max width size for presenter view
    # scv.settings.max_memory = args.ram  # for beautified visualization
    
    adata = read_loom(loom_path=args.input, clusterColumn=args.clusterColumn, stageColumn=args.stageColumn, batchColumn=args.batchColumn)
    plot_proportions(adata=adata, clusterColumn=args.clusterColumn)
    
    adata = preprocess_anndata(adata=adata)
    adata = calc_moments(adata=adata)
    adata = calc_velocity(adata=adata, velocityMode=args.velocityMode, ncores=args.ncores)
    plot_velocity(adata=adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
    

    # Calculate velocity pseudotime and paga graph and plot
    velocity_pseudotime(adata, dpi=args.dpi)
    adata, paga_df = paga(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
    
    
    # Run dynamical or deterministic models and plot genes
    if args.velocityMode == 'deterministic':
        if args.geneDiscovery == True:
            # Identify and plot genes that have cluster-specific differential velocity expression
            top_cluster_genes = ident_top_genes(adata, args.clusterColumn)
            plot_genes(adata=adata, genes_dict=top_cluster_genes, clusterColumn=args.clusterColumn, plot_dir = 'gene_discovery/', dpi=args.dpi)
            plot_differentiation(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
        # Plot velocity for manual GOI
        if args.genes is not None:
            manual_genes = {args.genes[i]: args.genes[i] for i in range(0, len(args.genes))}
            plot_genes(adata=adata, genes_dict=manual_genes, clusterColumn=args.clusterColumn, plot_dir = 'goi/', dpi=args.dpi)
    
    if args.velocityMode == 'dynamical':
        adata = latent_time(adata, dpi=args.dpi)
        if args.geneDiscovery == True:
            plot_latent_time_heatmap(adata, clusterColumn = args.clusterColumn, stageColumn = args.stageColumn, forceCol = args.forceCol, batchColumn = args.batchColumn, plot_dir='gene_discovery_dynamical/')
            # Identify and plot genes that have cluster-specific differential velocity expression
            top_cluster_genes = ident_top_genes_dynamical(adata, clusterColumn = args.clusterColumn)
            plot_genes_dynamical(adata, genes_dict=top_cluster_genes, clusterColumn=args.clusterColumn, plot_dir = 'gene_discovery_dynamical/', dpi=args.dpi)
        # Plot velocity for manual GOI
        if args.genes is not None:
            manual_genes = {args.genes[i]: args.genes[i] for i in range(0, len(args.genes))}
            plot_genes_dynamical(adata=adata, genes_dict=manual_genes, clusterColumn=args.clusterColumn, plot_dir = 'goi_dynamical/', dpi=args.dpi)
        
        print('Calculate differential kinetics and recompute velocity')
        adata = calc_diff_kinetics(adata=adata, clusterColumn=args.clusterColumn, plot_dir="")
        adata = calc_velocity(adata=adata, velocityMode=args.velocityMode, ncores=args.ncores, diffKinetics=True, groupby=args.clusterColumn)
        plot_velocity(adata=adata, clusterColumn=args.clusterColumn, plot_dir="diff_kinetics/", dpi=args.dpi)
    
    adata.write(args.output + '_scvelo.h5ad')
    adata.obs.to_csv(args.output + '_metadata.csv')

if __name__ == '__main__':
    sys.exit(main())


# # Generate test data
# args = parse_args(args)
# adata = read_loom(args.input, args.clusterColumn)
# adata = adata[adata.obs.index[0:1000], adata.var.index[0:5000]]
# adata.write_loom('test_loom.loom', write_obsm_varm=True)

# # set args for interactive testing
# args = ['-i', 'seurat_merged.loom', '-o', 'test.h5ad', '-m', 'dynamical', '-c', 'seurat_clusters', '-s', 'stage', '-b', 'run', '-n', '4', '--genes', 'NET1', 'ZC3HC1', '--differentialKinetics', 'True']
# main(args)