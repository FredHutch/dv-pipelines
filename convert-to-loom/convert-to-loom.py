#!/bin/python3
import sys
import os
import numpy as np
import pandas as pd
import scanpy.api as sc
import anndata
from anndata import AnnData

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

working_dir = os.getcwd()

adata = sc.read(os.path.join(working_dir, 'gene_count_head.txt'), cache=True).transpose()
print('adata variables:')
print(dir(data))
df_cell = pd.read_csv(os.path.join(working_dir, 'cell_annotate.csv'), delimiter='\t')
df_gene = pd.read_csv(os.path.join(working_dir, 'gene_annotate.csv'), delimiter="\t")
df_cell.index = df_cell["sample"]
df_gene.index = df_gene["gene_id"]
adata.obs = df_cell
adata.var = df_gene
# save the loom file
adata.write_loom("output.loom")
