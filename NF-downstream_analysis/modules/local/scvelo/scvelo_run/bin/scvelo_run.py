#!/opt/conda/envs/scvelo/bin/python

import os
import sys
import argparse
import scvelo as scv
import pdfkit


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    parser.add_argument('-m', '--velocityMode', help="Method for calculating RNA velocity. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.", default='dynamical', metavar='')
    parser.add_argument('-x', '--clusterColumn', help="Name of cluster column.", default='Clusters', metavar='')
    parser.add_argument('-g', '--genes', help="Genes of interest to plot on velocity.", nargs='+')
    parser.add_argument('-c', '--ncores', help="Number of cores used for parallelisation.", metavar='')
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
    print(ncores)
    if velocityMode not in ['dynamical', 'deterministic', 'stochastic']:
        Exception(f"'--velocityMode': '{velocityMode}' is not valid. Must be set to either: 'dynamical', 'deterministic', or 'stochastic'.")

    if velocityMode == 'dynamical':
        scv.tl.recover_dynamics(adata, n_jobs=ncores)

    scv.tl.velocity(adata, mode=velocityMode)
    scv.tl.velocity_graph(adata)
    return(adata)

def plot_velocity(adata, clusterColumn):
    scv.pl.velocity_embedding(adata, color=clusterColumn, arrow_length=3, arrow_size=2, dpi=120, save='embedding.png')
    scv.pl.velocity_embedding_grid(adata, color=clusterColumn, basis='umap', save='embedding_grid.png')
    scv.pl.velocity_embedding_stream(adata, color=clusterColumn, basis='umap', save='embedding_stream.png')


# identify genes which explain the vector field in a given lineage
def ident_genes(adata, clusterColumn, min_corr=.3):
    scv.tl.rank_velocity_genes(adata, groupby=clusterColumn, min_corr=min_corr)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    return(df)

def plot_genes(adata, args):
    if args.genes is not None:
        scv.pl.velocity(adata, args.genes, ncols=2)
#     df = ident_genes(adata, args.clusterColumn)
#     for column in df:
#         scv.pl.scatter(adata, df[column][:5], ylabel=column)
#         scv.pl.velocity(adata, df[column][:5], ncols=5, ylabel=column, dpi=240)

def plot_differentiation(adata, clusterColumn):
    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save='differentiation_scatter.png')
    
    df = adata.obs.groupby(clusterColumn)[keys].mean().T
    styled_table = df.style.background_gradient(cmap='coolwarm', axis=1)
    html = styled_table.render()
    pdfkit.from_string(html, 'figures/cluster_differentiation.png')



def main(args=None):
    print(args.ncores)
    args = parse_args(args)
    # check_args(args)
    
    adata = read_loom(args.input, args.clusterColumn)
    
    plot_proportions(adata, args.clusterColumn)
    preprocess_anndata(adata)
    calc_moments(adata)
    calc_velocity(adata, args.velocityMode, args.ncores)
    plot_velocity(adata, args.clusterColumn)
    
    plot_genes(adata, args)
    plot_differentiation(adata, args.clusterColumn)
    
    latent_time(adata)
        
    return(args, adata)


if __name__ == "__main__":
    sys.exit(main())



# set args for interactive testing
# args = ['-i', '../output/NF-downstream_analysis_stacas/scvelo/seurat_intersect_loom/seurat_merged.loom', '-vm', 'deterministic']
# args = ['-i', 'test_loom.loom', '-vm', 'deterministic']

# args, adata = main(args)

# # Generate test data
# args = parse_args(args)
# adata = read_loom(args.input)
# adata = adata[adata.obs.index[0:1000], adata.var.index[0:5000]]
# adata.write_loom('test_loom.loom', write_obsm_varm=True)

