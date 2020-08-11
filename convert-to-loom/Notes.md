# Pulling the Datasets


## Repro on a Laptop

1. Create a new env for python 3.8 in conda/virtualenv

```bash
conda create -n converttoloom38 python=3.8 anaconda 
```

2. Install the the required packages using pip
```bash
conda activate converttoloom38
python --version

python -m pip install scanpy pandas

python -m pip install cd /Users/dnambi/Documents/GitHub/data-vizualization-center/hdf5/pubweb

jupyter notebook

conda converttoloom38

```


What's going on with the array?
adata.X is an numpy.ndarray

cell_annotate.sample = `sci3-me-001.GTCGGAGTTTGAGGTAGAA`
gene_annotate.gene_id = `ENSMUSG00000051951.5`

Use awk to figure out which column maps to which. Whaaat? 

```python
adata.obs = df_cell
adata.var = df_gene
```

[Anndata](https://anndata.readthedocs.io/en/latest/concatenation.html#annotating-data-source-label-keys-and-index-unique) and 





Set up an Anaconda environment

```bash

conda activate 



```



## Copy in the Code

```bash

cd /Users/dnambi/Documents/GitHub/data-vizualization-center/hdf5
tar -czf "pubweb.tar.gz" pubweb/*

sshh

scp /Users/dnambi/Documents/GitHub/data-vizualization-center/hdf5/pubweb.tar.gz dnambi@rhino:/fh/fast/_HDC/user/dnambi/anndata
```

On Rhino:

```bash
cd /fh/fast/_HDC/user/dnambi/anndata
tar -xzf pubweb.tar.gz

rm pubweb.tar.gz
```

## This dataset

Use ```/fh/fast/_HDC/user/dnambi``` for storage for now

Use ```/fh/fast/_HDC/user/dnambi/anndata/cao1025``` for this sample script

Make the slurm command:


```bash
cd /fh/fast/_SR/dataviz/anndata/cao1025

# ask for lots of CPU as a proxy for a lot of memory
sbatch -c 8 -J convert-to-loom wrapper.sh

```
Get the anndata.py file:

```bash
scp -r dnambi@rhino:/fh/fast/_SR/dataviz/anndata/cao1025/anndata.py /Users/dnambi/Documents

scp -r dnambi@rhino:/fh/fast/_SR/dataviz/anndata/cao1025/gene_count_head.txt /Users/dnambi/Documents

```

Gene count has 1,251,379,870 rows. Let's try with 1 million rows.

```bash
head -n 1000000 gene_count.txt > gene_count_head.txt
```

Fix up the header:

```bash
sed -i 's/general//' gene_count.txt

sed -i 's/ general//' work

```


Try getting the numbers to align:
```bash
wc -l gene_annotate.csv #26184
wc -l cell_annotate.csv #2058653
```

count * count = 53,903,770,152, 53 billion, not. Unless we're dealing with a sparse matrix? 


Running the thing end to end:

```bash
cd 
sed -i 's/ general//g' gene_count.txt

sbatch -c 8 -J convert-to-loom-v3 wrapper.sh

```


Error:

```python
Traceback (most recent call last):
  File "convert-to-loom-3.py", line 22, in <module>
    adata.obs = df_cell
  File "/app/software/scanpy/1.4.6-foss-2019b-Python-3.7.4/lib/python3.7/site-packages/anndata/_core/anndata.py", line 801, in obs
    self._set_dim_df(value, "obs")
  File "/app/software/scanpy/1.4.6-foss-2019b-Python-3.7.4/lib/python3.7/site-packages/anndata/_core/anndata.py", line 750, in _set_dim_df
    value_idx = self._prep_dim_index(value.index, attr)
  File "/app/software/scanpy/1.4.6-foss-2019b-Python-3.7.4/lib/python3.7/site-packages/anndata/_core/anndata.py", line 763, in _prep_dim_index
    f"Length of passed value for {attr}_names is {len(value)}, but this AnnData has shape: {self.shape}"
ValueError: Length of passed value for obs_names is 2058652, but this AnnData has shape: (3, 1251379869)
```




## Other Sets

https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing

http://atlas.gs.washington.edu/mouse-atac/

http://atlas.gs.washington.edu/fly-atac/

http://atlas.gs.washington.edu/worm-rna/