import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import anndata as ad

# Plotting libraries
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Inline plotting settings
%matplotlib inline
sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

# Reset random seed
np.random.seed(5)

# Load expression data
expr = pd.read_csv('expr_hvg.csv', delimiter=',', index_col=0)
adata = ad.AnnData(X=expr)

# Perform PCA
sc.pp.pca(adata, n_comps=10)
pca_projections = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections)

# Determine multiscale space
ms_data = palantir.utils.determine_multiscale_space(dm_res)

# Compute neighbors and draw graph
sc.pp.neighbors(adata, n_neighbors=50)
sc.tl.draw_graph(adata, layout='fa')

# Plot the force-directed layout
sc.pl.draw_graph(adata, layout='fa')

# Extract force-directed layout coordinates
fa = pd.DataFrame(adata.obsm['X_draw_graph_fa'], index=adata.obs_names)
fa.columns = ['x', 'y']

# Define the starting cell (cell with the highest expression of photosynthesis gene)
start_cell = 'TACACGACAGACAAAT-1'

# Run Palantir
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500)

# Plot Palantir results on the force-directed layout
palantir.plot.plot_palantir_results(pr_res, fa)