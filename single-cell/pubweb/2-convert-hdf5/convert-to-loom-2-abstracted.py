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
import argparse
import json

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

logging.basicConfig(filename='convert-to-loom.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Convert input to HDF5")
parser.add_argument('--params')
parser.add_argument('--matrix')
parser.add_argument('--obs')
parser.add_argument('--var')
parser.add_argument('--output')

args = parser.parse_args()

with open(args.params, 'r') as f:
    json_params = json.loads(f.read())

try:
    working_dir = os.getcwd()
    # get the parameter object
    params = json_params.get('parameters').get('input').get('anndata-spec')
    # get the matrix file
    logging.info('Now reading the matrix market file')
    adata = sc.read_mtx(os.path.join(working_dir, args.matrix)).transpose()
    logging.debug('adata variables:')
    logging.info(dir(adata))
    # get the var dataframe
    logging.info('Now reading var dataframe')
    df_var = pd.read_csv(os.path.join(working_dir, args.var),
                         delimiter=params.get('var').get('delimiter'))
    logging.debug('var columns are {}'.format(list(df_var.columns)))
    logging.debug('Now mapping pandas indexes')
    df_var.index = df_var[params.get('var').get('index_col')]
    df_var.index.name = None
    logging.debug('Now mapping var')
    adata.var = df_var
    # get the obs data
    logging.info('Now reading obs dataframe')
    df_obs = pd.read_csv(os.path.join(working_dir, args.obs),
                         delimiter=params.get('obs').get('delimiter'))
    logging.debug('obs columns are: {}'.format(list(df_obs.columns)))
    df_obs.index = df_obs[params.get('obs').get('index_col')]
    df_obs.index.name = None
    logging.debug('Now mapping obs')
    adata.obs = df_obs
    # map in observational matrices
    logging.debug('Now mapping the obsm')
    adata.obsm = None
    for k in params.get('obsm').keys():
        adata.obsm[k] = adata.obs[params.get('obsm').get(k)].values
        logging.debug(f"obsm[{k}] = obs[{params.get('obsm').get(k)}]")
    # now drop the index columns, they're already in the index
    df_var.drop(params.get('var').get('index_col'), axis=1, inplace=True)
    df_obs.drop(params.get('obs').get('index_col'), axis=1, inplace=True)
    # save the hdf5 file
    logging.info(f'Now writing hdf5 file to {args.output}')
    adata.write(args.output)
except Exception as E:
    tracer = traceback.format_exc()
    type_e, value, ignore = sys.exc_info()
    logging.error(f"ERROR: {repr(E)} {type_e} value {tracer}")
