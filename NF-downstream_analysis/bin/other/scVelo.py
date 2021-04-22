import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
%load_ext rpy2.ipython

sample_one = anndata.read_loom("sample_one.loom")
....
sample_n = anndata.read_loom("sample_n.loom")

sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters_obs.csv")


sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["x"])]