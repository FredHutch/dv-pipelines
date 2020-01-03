#!/usr/bin/env nextflow

INPUT_CH = Channel.from([[
  file(params.input.exp), 
  file(params.input.col), 
  file(params.input.row), 
  UUID.randomUUID().toString().substring(0,7)
]])

INPUT_CH.into {
  TENX_INPUT_CH
  CDS3_INPUT_CH
}

process TENX_LOAD_PR {

  publishDir "$params.output.folder/format/tenx"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file('exp.rds'), file('col.rds'), file('row.rds'), val(ppid) from TENX_INPUT_CH

  output:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv') into TENX_LOAD_CH

  """
  #!/usr/bin/env Rscript 
  library(Matrix)
  library(Seurat)
  expression_matrix <- readRDS("exp.rds")
  cell_metadata <- readRDS("col.rds")
  gene_annotation <- readRDS("row.rds")
  writeMM(expression_matrix, 'matrix.mtx')
  write.table(row.names(cell_metadata), "barcodes.tsv", quote=FALSE, sep='\t', row.names=F, col.names=F)
  write.table(cell_metadata, "cells.tsv", quote=FALSE, sep='\t', row.names=T, col.names=T)
  write.table(gene_annotation, "genes.tsv", quote=FALSE, sep='\t', row.names=T, col.names=F)
  """
}

TENX_LOAD_CH.into {
  SEURAT_INPUT_CH 
  SCANPY_INPUT_CH
}

process CDS3 {

  publishDir "$params.output.folder/format/cds3"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file('exp.rds'), file('col.rds'), file('row.rds'), val(ppid) from CDS3_INPUT_CH

  output:
    file "data.rds"

  """
  monocle3 create -F cds3 --cell-metadata=col.rds --gene-annotation=row.rds --expression-matrix=exp.rds data.rds
  """
}

process SCANPY_PR {

  publishDir "$params.output.folder/format/h5ad"
  container "quay.io/biocontainers/scanpy-scripts:0.0.3--py37_1"

  input:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv') from SCANPY_INPUT_CH

  output:
    file "data.h5ad"

  """
  scanpy-read-10x.py -d \$PWD/ -o data.h5ad -F anndata
  """
}

process SEURAT_PR {

  publishDir "$params.output.folder/format/seurat/v2"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv') from SEURAT_INPUT_CH

  output:
    file "data.rds" into SEURAT_LOAD_CH 

  """
  Rscript /usr/local/bin/seurat-read-10x.R -d \$PWD -o seurat-tmp.rds
  Rscript /usr/local/bin/seurat-create-seurat-object.R -i seurat-tmp.rds -o data.rds
  """
}

process SEURAT_UPDATE_PR { 

  publishDir "$params.output.folder/format/seurat/v3"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    file("data.rds") from SEURAT_LOAD_CH

  output:
    file("data.rds") into SEURAT_UPDATE_CH

  """
  echo 'library(Seurat)' >> tmp.R
  echo 'data <- readRDS("data.rds")' >> tmp.R
  echo 'data <- UpdateSeuratObject(data)' >> tmp.R
  echo 'saveRDS(data, "data.rds")' >> tmp.R
  R -f tmp.R
  """
}

// SEURAT_UPDATE_CH.into{
//   SCE_INPUT_CH
//   LOOM_INPUT_CH
// }

// process SCE_PR {
   
//   publishDir "$params.output.folder/format/sce"
//   container "quay.io/biocontainers/r-seurat:3.0.2--r36h0357c0b_1"

//   input:
//     file('data.rds') from SEURAT_UPDATE_CH

//   output:
//     file "data.rds" 

//   """
//   #!/usr/bin/env Rscript 
//   library(Seurat)
//   seuratObj <- readRDS("data.rds")
//   sce <- as.SingleCellExperiment(seuratObj)
//   saveRDS(sce, "data.rds")
//   """
// }

// process LOOM_PR {

//   publishDir "$params.output.folder/format/loom"
//   container "quay.io/biocontainers/r-seurat:3.0.2--r36h0357c0b_1"

//   input:
//     file('data.rds') from LOOM_INPUT_CH

//   output:
//     file "data.rds" 

//   """
//   #!/usr/bin/env Rscript 
//   library(Seurat)
//   seuratObj <- readRDS("data.rds")
//   loom <- as.loom(seuratObj)
//   saveRDS(loom, "data.loom")
//   """
// }
