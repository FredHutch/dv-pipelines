#!/bin/python3
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from anndata import AnnData
import logging
import traceback

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

logging.basicConfig(filename='convert-to-loom.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

try:
    working_dir = os.getcwd()
    logging.info('Now reading the matrix market file')
    adata = sc.read_mtx(os.path.join(working_dir, 'gene_count.txt')).transpose()
    # adata = sc.read(os.path.join(working_dir, 'gene_count_head.txt'), cache=True).transpose()
    logging.debug('adata variables:')
    logging.info(dir(adata))
    logging.info('Now reading cell_annotate.csv')
    df_cell = pd.read_csv(os.path.join(working_dir, 'cell_annotate.csv'), delimiter=',')
    logging.info('cell_annotate columns are: {}'.format(list(df_cell.columns)))
    logging.info('Now reading gene_annotate.csv')
    df_gene = pd.read_csv(os.path.join(working_dir, 'gene_annotate.csv'), delimiter=',')
    logging.info('gene_annotate columns are {}'.format(list(df_gene.columns)))
    logging.info('Now mapping pandas indexes')
    df_cell.index = df_cell["sample"]
    df_gene.index = df_gene["gene_id"]
    df_cell.index.name = None
    df_gene.index.name = None
    logging.info('Now mapping var')
    adata.var = df_gene
    logging.info('Now mapping obs')
    adata.obs = df_cell # based off of index
    logging.info('Now mapping the obsm')
    adata.obsm = None
    adata.obsm['Tsne_main_cluster'] = adata.obs[['Main_cluster_tsne_1', 'Main_cluster_tsne_2']].values
    adata.obsm['Tsne_sub_cluster'] = adata.obs[['Sub_cluster_tsne_1', 'Sub_cluster_tsne_2']].values
    adata.obsm['Umap_main_trajectory'] = adata.obs[['Main_trajectory_umap_1', 'Main_trajectory_umap_2', 'Main_trajectory_umap_3']].values
    adata.obsm['Umap_main_trajectory_refined'] = adata.obs[['Main_trajectory_refined_umap_1', 'Main_trajectory_refined_umap_2', 'Main_trajectory_refined_umap_3']].values
    adata.obsm['Umap_sub_trajectory'] = adata.obs[['Sub_trajectory_umap_1', 'Sub_trajectory_umap_2']].values
    # now drop the index columns, they're already in the index
    df_gene.drop('gene_id', axis=1, inplace=True)
    df_cell.drop('sample', axis=1, inplace=True)
    # save the loom file
    logging.info('Now writing hdf5 file using default')
    adata.write('output.hdf5')
except Exception as E:
    tracer = traceback.format_exc()
    type_e, value, ignore = sys.exc_info()
    logging.error(f"ERROR: {repr(E)} {type_e} value {tracer}")
