#!/bin/python3
import sys
import os
import pandas as pd
import scanpy as sc
import anndata
from anndata import AnnData

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

working_dir = os.getcwd()

adata = sc.read(os.path.join(working_dir, 'gene_count.txt'),
                cache=True,
                delimiter=' ').transpose()
print('adata variables:')
print(dir(adata))
df_cell = pd.read_csv(os.path.join(working_dir, 'cell_annotate.csv'), delimiter=',')
df_gene = pd.read_csv(os.path.join(working_dir, 'gene_annotate.csv'), delimiter=',')
df_cell.index = df_cell["sample"]
df_gene.index = df_gene["gene_id"]
adata.obs = df_cell
adata.var = df_gene
# save the loom file
adata.write_loom("output.loom")