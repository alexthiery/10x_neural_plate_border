#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import scvelo as scv


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    parser.add_argument('-m', '--velocityMode', help="Method for calculating RNA velocity. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.", default='dynamical', metavar='')
    parser.add_argument('-x', '--clusterColumn', help="Name of cluster column.", default='Clusters', metavar='')
    parser.add_argument('-g', '--genes', help="Genes of interest to plot on velocity.", nargs='+')
    parser.add_argument('-c', '--ncores', help="Number of cores used for parallelisation.", metavar='')
    parser.add_argument('-d', '--dpi', type=int, help='Set DPI for plots.', default='240')
    return parser.parse_args(args)
    

# def check_args(args=None):
#     if os.path.isdir(args.input):
#         if args.loomOutput is None:
#             raise Exception("'--loomOutput' must be specified when '--input' is a directory containing loom files to concatenate.")
#     elif not os.path.isfile(args.input):
#         raise Exception(f"'--input': '{args.input}' is not a valid path.")


# Read in loom data
def read_loom(loom_path, clusterColumn):
    adata = scv.read(loom_path)
    # set cluster column as categorical variable
    adata.obs[clusterColumn]=adata.obs[clusterColumn].astype('category')
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
def calc_velocity(adata, velocityMode, ncores):
    if velocityMode not in ['dynamical', 'deterministic', 'stochastic']:
        Exception(f"'--velocityMode': '{velocityMode}' is not valid. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.")

    if velocityMode == 'dynamical':
        scv.tl.recover_dynamics(adata, n_jobs=ncores)

    scv.tl.velocity(adata, mode=velocityMode)
    scv.tl.velocity_graph(adata)
    return(adata)

def plot_velocity(adata, clusterColumn, threshold=.1, dpi=240):
    scv.pl.velocity_embedding(adata, color=clusterColumn, arrow_length=3, arrow_size=2, save='embedding.png', dpi=dpi)
    scv.pl.velocity_embedding_grid(adata, color=clusterColumn, basis='umap', save='embedding_grid.png', dpi=dpi)
    scv.pl.velocity_embedding_stream(adata, color=clusterColumn, basis='umap', save='embedding_stream.png', dpi=dpi)
    scv.pl.velocity_graph(adata, threshold=threshold, color=clusterColumn, basis='umap', save='graph.png', dpi=dpi)

# plot expression of genes of interest across velocity UMAPs
def plot_genes(adata, genes_dict, dpi=240, prefix=""): # genes_dict is a key value pair with {name:gene(s)}
    for key, value in genes_dict.items():
        scv.pl.velocity(adata, value, save=prefix+key+'.png', dpi=dpi)

# identify genes which explain the vector field in a given lineage
def ident_top_genes(adata, clusterColumn, min_corr=.3, n_genes=5):
    scv.tl.rank_velocity_genes(adata, groupby=clusterColumn, min_corr=min_corr)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    
    top_cluster_genes = {'Cluster-'+column: df[column][:n_genes] for column in df}
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

def latent_time(adata, clusterColumn, dpi=240):
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png', dpi=dpi)
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color=clusterColumn, n_convolve=100, save='heatmap.png', dpi=dpi)
    return(adata)

def main(args=None):
    
    args = parse_args(args)
    
    # check_args(args)
    
    adata = read_loom(loom_path=args.input, clusterColumn=args.clusterColumn)
    
    plot_proportions(adata=adata, clusterColumn=args.clusterColumn)
    
    adata = preprocess_anndata(adata=adata)
    
    adata = calc_moments(adata=adata)
    
    adata = calc_velocity(adata=adata, velocityMode=args.velocityMode, ncores=args.ncores)
    plot_velocity(adata=adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
    
    # Plot velocity for manual GOI
    if args.genes is not None:
        manual_genes = {args.genes[i]: args.genes[i] for i in range(0, len(args.genes))}
        plot_genes(adata=adata, genes_dict=manual_genes, dpi=args.dpi)
        
    # Identify and plot genes that have cluster-specific differential velocity expression
    top_cluster_genes = ident_top_genes(adata, args.clusterColumn)
    plot_genes(adata=adata, genes_dict=top_cluster_genes, dpi=args.dpi)
    
    plot_differentiation(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
    
    Calculate velocity pseudotime and plot
    velocity_pseudotime(adata, dpi=args.dpi)
    
    adata, paga_df = paga(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
    adata = latent_time(adata, clusterColumn=args.clusterColumn, dpi=args.dpi)
        
    return(args, adata)

if __name__ == "__main__":
    sys.exit(main())


# # Generate test data
# args = parse_args(args)
# adata = read_loom(args.input, args.clusterColumn)
# adata = adata[adata.obs.index[0:1000], adata.var.index[0:5000]]
# adata.write_loom('test_loom.loom', write_obsm_varm=True)

# set args for interactive testing
# args = ['-i', '../output/NF-downstream_analysis_stacas/scvelo/seurat_intersect_loom/seurat_merged.loom', '-m', 'deterministic', '-x', 'Clusters']
# args = ['-i', 'test_loom.loom', '-m', 'deterministic', '-x', 'Clusters', '-c', '2', '-g', 'BCLAF1', 'AHI1', '--dpi', '80']

# args, adata = main(args)
