#! /usr/bin/env python

import argparse
from pubweb.singlecell import CellRanger

parser = argparse.ArgumentParser(description="Invokes Cellranger HDF5 convert")
parser.add_argument('--name')
parser.add_argument('--input')
parser.add_argument('--output')
parser.add_argument('--species')

args = parser.parse_args()

CellRanger(
    inputFolder=args.input,
    outputFolder=args.output,
    datasetName=args.name,
    species=args.species)
