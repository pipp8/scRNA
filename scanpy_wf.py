#! /usr/bin/python3

import string
import numpy as np
import sys
from pathlib import Path
from  datetime import  datetime

# import library
import pandas as pd
import scanpy as sc

import anndata
import openpyxl
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score




def runTest(input: string):
    inputPath = Path(input)
    ext = inputPath.suffix

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    #sc.logging.print_header()
    #sc.settings.set_figure_params(dpi=80, facecolor='white')

    # save time usage ####
    time_sc = pd.DataFrame(index=["loading", "find_mit_gene", "filter", "normalization", "hvg",
                                  "scaling", "PCA", "t-sne", "umap", "louvain", "leiden"],
                           columns=["time_sec"])

    # data ####
    start_time = datetime.now()

    if (ext == '.h5ad'):
        adata = sc.read_h5ad(input)
    elif (ext == '.h5'):
        adata = sc.read_10x_h5(input)
    else:
        print(f'This script can read only h5ad | TenX.h5 files')
        exit(-2)


    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    print( "Input size: ", adata.shape, flush=True)
    adata

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Loading data Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[0, 0] = time_elapsed

    start_time = datetime.now()
    # Show those genes that yield the highest fraction of counts in each single cell, across all cells.
    sc.pl.highest_expr_genes(adata, n_top=20, )
    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Selecting top 20 genes Time Elapsed:", time_elapsed, flush=True)

    # find mitocondrial genes ####
    start_time = datetime.now()
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Finding MT genes Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[1, 0] = time_elapsed

    # filter data sulle cellule, profondità di sequenziamento####
    start_time = datetime.now()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    adata
    print(f"{end_time}: Filtering Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[2, 0] = time_elapsed

    # normalization ####
    start_time = datetime.now()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Normalization Time Elapsed:", time_elapsed)
    time_sc.iloc[3, 0] = time_elapsed

    # Identification of highly variable features (feature selection) ####
    start_time = datetime.now()

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,
                                n_top_genes = 1000)
    sc.pl.highly_variable_genes(adata)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    adata
    print(f"{end_time}: selecting higly variable genes Time Elapsed:", time_elapsed)
    time_sc.iloc[4, 0] = time_elapsed

    df = pd.DataFrame(adata.var.highly_variable, columns=['hvg'])
    df.to_excel(inputPath.with_name( inputPath.suffix + '_hvg.xlsx'), index=False)



    # Scaling the data ####
    start_time = datetime.now()

    sc.pp.scale(adata, max_value=10)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Scaling Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[5, 0] = time_elapsed

    print( f"{end_time}: PCA input size: ", adata.shape, flush=True)
    adata.write_csvs(inputPath.with_name( inputPath.stem + '_scaled'), False)

    # PCA ####
    start_time = datetime.now()

    sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca(adata, color='CST3')
    sc.pl.pca_variance_ratio(adata, log=True)
    # adata.write(results_file) capire come definirlo
    adata

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: PCA Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[6, 0] = time_elapsed
    adata.write_csvs(inputPath.with_name( inputPath.stem + '_PCA'), False)

    # t-sne ####
    start_time = datetime.now()

    sc.tl.tsne(adata)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: T-SNE Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[7, 0] = time_elapsed

    # UMAP ####
    start_time = datetime.now()
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: UMAP Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[8, 0] = time_elapsed

    # louvain ####
    start_time = datetime.now()

    sc.tl.louvain(adata, resolution = 0.13)

    end_time = datetime.now()
    time_elapsed = end_time - start_time
    print(f"{end_time}: Louvain Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[9, 0] = time_elapsed

    true_labels = adata.obs['Sample'].astype(str)
    predicted_labels = adata.obs['louvain'].astype(str)



    # Compute the ARI
    ari_score = adjusted_rand_score(true_labels, predicted_labels)

    # Print the ARI score
    print(f"{end_time}: Adjusted Rand Index (ARI):", ari_score)

    cluster_labels = adata.obs['louvain']

    # Calculate silhouette scores
    silhouette_avg = silhouette_score(adata.X, cluster_labels)

    # Print the average silhouette score
    print(f"{end_time}: Average Silhouette Score:", silhouette_avg)


    # leiden ####
    start_time = datetime.now()

    sc.tl.leiden(adata, resolution = 0.13)

    end_time = datetime.now() 
    time_elapsed = end_time - start_time
    print(f"{end_time}: Leiden Time Elapsed:", time_elapsed, flush=True)
    time_sc.iloc[10, 0] = time_elapsed

    true_labels = adata.obs['Sample'].astype(str)
    predicted_labels = adata.obs['leiden'].astype(str)

    # Compute the ARI
    ari_score = adjusted_rand_score(true_labels, predicted_labels)

    # Print the ARI score
    print(f"{end_time}: Adjusted Rand Index (ARI):", ari_score)

    cluster_labels = adata.obs['leiden']

    # Calculate silhouette scores
    silhouette_avg = silhouette_score(adata.X, cluster_labels)

    # Print the average silhouette score
    print(f"{end_time}: Average Silhouette Score:", silhouette_avg)

    time_sc

    adata.write_csvs(inputPath.with_name( inputPath.suffix + '_scanpy'), False)





if __name__ == "__main__":
    argNum = len(sys.argv)
    if (argNum != 2):
        print( f"{Path(sys.argv[0]).name} inputData.h5ad")
        exit(-1)
    else:
        runTest(sys.argv[1])
