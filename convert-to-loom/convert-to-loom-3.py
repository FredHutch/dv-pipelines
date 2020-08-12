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
adata.obs = df_cell # based off of index
adata.var = df_gene
# save the loom file
adata.write_loom("output.loom")


# gene count - cell by gene, based on index. They're exported using the same index
# annotations on the columns
# gene index annotations on the rows
# other thing to look at - anndata. CDS objects (cell dataset, in R)
# h5ad, loom, cds, and Seurat. And now our pubweb format.


top row - dimensions of the dataset. It's a sparse matrix.
It's mostly zeroes. 
58 - 1 1 -> 
# matrix market format, https://math.nist.gov/MatrixMarket/formats.html
First column is row, second column is column, the third is value. Index is 1 based
# wilfred is struggling to convert this to a dense matrix. -> they're straight up data tables.

HDF5 -either sparse matrices or strongly typed big-ass lists. 



