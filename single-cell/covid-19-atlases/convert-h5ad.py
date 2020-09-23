#! /usr/bin/env python
import argparse
import anndata
import h5py
import os


parser = argparse.ArgumentParser(description="Converts anndata h5ad to hdf5")
parser.add_argument('--input')
parser.add_argument('--output')

args = parser.parse_args()

print(f"Converting {args.input} to {args.output}")
adata = anndata.read_h5ad(args.input, backed='r')
adata.write(args.output)

print(f"Changing index groups for {args.output}")
with h5py.File(args.output, "r+") as data:
    if 'index' in data['var'].keys():
        print('Moving /var/index group')
        data.move('var/index', 'var/_index')
    if 'index' in data['obs'].keys():
        print('Moving /obs/index group')
        data.move('obs/index', 'obs/_index')
