import scvelo as scv
import numpy as np
import os
import anndata as ad
import pandas as pd
import re

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")

def load():
    # Load the loom file
    adata = scv.read("/home/cwlee/project/DGIST_REC/scvelo/arabidopsis.loom", cache=True)
    adata.var_names_make_unique()

    # Add Seurat cluster information to adata
    sample_obs = pd.read_csv("/home/cwlee/project/DGIST_REC/scvelo/final2/cellID_obs.csv")
    arabidopsis_IND = list(map(lambda x: re.sub("(.*)-.*","counts:\\1x",x) ,sample_obs.iloc[:,0].values))
    sample_obs.loc[:,"loomstyle"] = arabidopsis_IND
    adata = adata[np.isin(adata.obs.index,sample_obs["loomstyle"])]

    # Load cluster information
    stages = pd.read_csv('/home/cwlee/project/DGIST_REC/scvelo/final2/clusters.csv')
    arabidopsis_IND = list(map(lambda x: re.sub("(.*)-.*","counts:\\1x",x) ,stages.iloc[:,0].values))
    stages.loc[:,"loomstyle"] = arabidopsis_IND
    stages.index = stages["loomstyle"]
    adata.obs['Clusters'] = stages.loc[adata.obs['Clusters'].index.values,:]['x']

    # Load UMAP embeddings
    fa_coordinate = pd.read_csv("/home/cwlee/project/DGIST_REC/scvelo/final2/cell_embeddings.csv")
    arabidopsis_IND = list(map(lambda x: re.sub("(.*)-.*","counts:\\1x",x) ,fa_coordinate.iloc[:,0].values))
    fa_coordinate.loc[:,"loomstyle"] = arabidopsis_IND
    fa_coordinate.index = fa_coordinate["loomstyle"]
    adata.obsm['X_fa'] = np.array([list(fa_coordinate.loc[adata.obs['Clusters'].index.values,:]['palantir_1'].values),
                                     list(fa_coordinate.loc[adata.obs['Clusters'].index.values,:]['palantir_2'].values)]).transpose()
    return adata

adata = load()

# Preprocessing
scv.pp.filter_and_normalize(adata, min_shared_counts=5, n_top_genes=500)
scv.pp.moments(adata, n_pcs=15, n_neighbors=30)

# Velocity analysis
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Set cluster colors
adata.uns['Clusters_colors'] = ["#2B7A0B","#7DCE13","#EAE509"]

# Plot velocity embedding stream
scv.pl.velocity_embedding_stream(adata, basis='fa', color=['Clusters'], legend_loc="right margin", save='./scvelo/final2/scvelo.png')