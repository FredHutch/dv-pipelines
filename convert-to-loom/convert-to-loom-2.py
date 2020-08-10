#!/bin/python3
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

working_dir = os.getcwd()

adata = sc.read(os.path.join(working_dir, 'gene_count_head.txt'),
                cache=True,
                delimiter=' ').transpose()
print('adata variables:')
print(dir(adata))