#!/bin/bash
source /app/lmod/lmod/init/profile

# get the data
ml Python/3.8.2-GCCcore-9.3.0
ml scanpy


python -m pip install /fh/fast/_SR/dataviz/anndata/pubweb

echo "Now replacing the unneeded column header"
sed -i 's/ general//' gene_count.txt

echo "Now calling the Python script"
# call the python script
python convert-to-loom-3.py
