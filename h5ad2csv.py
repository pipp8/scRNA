#! /usr/bin/python3

import string
import numpy as np
import sys
from pathlib import Path

# import library
import pandas as pd
import scanpy as sc
import time
import anndata
import openpyxl
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score




def runTest(input: string):
    inputPath = Path(input)

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    #sc.logging.print_header()
    #sc.settings.set_figure_params(dpi=80, facecolor='white')


    # data ####

    adata = sc.read_h5ad(input)
    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

    adata = adata.transpose()
    
    print( "input size: ", adata.shape, flush=True)
    
    adata.write_csvs(inputPath.with_name( inputPath.stem + '_CSV'), False)






if __name__ == "__main__":
    argNum = len(sys.argv)
    if (argNum != 2):
        print( f"{Path(sys.argv[0]).name} inputData.h5ad")
        exit(-1)
    else:
        if (Path(sys.argv[1]).suffix != '.h5ad'):
            print(f'This script can read only h5ad input files')
            exit(-2)
        else:
            runTest(sys.argv[1])
