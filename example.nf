#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process construct {
  // S3 Path
}

process run {
    container "quay.io/biocontainers/scanpy"
    cpus 1
    memory '512 MB'

  input:
    path inputLoomFile

  output:
    file 'data.loom'

  script:

  """
  #!/usr/bin/python
  import scanpy as sc
  adata = sc.read_loom($loomfile)
  adata.var_names_make_unique()
  sc.pp.recipe_seurat(adata)
  sc.pp.pca(adata, copy=True)
  sc.write("data.h5", adata)
  """
}

process destruct {
  
}

workflow {
  construct()
  run()
  destruct()
}