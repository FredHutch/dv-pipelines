#!/usr/bin/env nextflow

process main {

"""
import scanpy as sc
adata = sc.read_loom("/Users/zager/Desktop/HPSI_human_cerebral_organoids_homo_sapiens.loom")
adata.var_names_make_unique()
sc.pp.recipe_seurat(adata)
sc.pp.pca(adata, copy=True)
sc.write("/Users/zager/Desktop/data.h5", adata)
"""

}