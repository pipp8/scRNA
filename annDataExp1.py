import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


print(ad.__version__)

counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)

# AnnData object with n_obs × n_vars = 100 (cells) × 2000 (genes)
adata
# <Compressed Sparse Row sparse matrix of dtype 'float32'
# 	with 126458 stored elements and shape (100, 2000)>
adata.X

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])

adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
