#!/usr/bin/env nextflow

INPUT_CH = Channel.from([[
  file(params.output.folder+ '/format/h5ad/data.h5ad'), 
  UUID.randomUUID().toString().substring(0,7)
]])

process SCANPY_NORM_PR { 
  publishDir "$params.output.folder/analysis/scanpy"
  container "zager/scanpy:3"

  input:
    set file("data.h5ad"), val(ppid) from INPUT_CH
    each scale_factor from params.scanpy.normalise.scale_factor

  output:
    set file("${ppid + '-' + pid}.h5ad"), file("${ppid +'-'+ pid}.txt"), val(pid) into SCANPY_NORM_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  scanpy-cli norm -t ${scale_factor} -f anndata data.h5ad ${ppid + '-' + pid}.h5ad 
  echo "scanpy filter-cells" >> ${ppid +'-'+ pid}.txt
  """
}

process SCANPY_PCA_PR { 
  publishDir "$params.output.folder/analysis/scanpy"
  container "zager/scanpy:3"

  input:
    set file("data.h5ad"), file("ledger.txt"), val(ppid) from SCANPY_NORM_CH
    each n_comps from params.scanpy.pca.n_comps
    each svd_solver from params.scanpy.pca.svd_solver
    each use_all from params.scanpy.pca.use_all

  output:
    set file("${ppid + '-' + pid}.h5ad"), file("${ppid +'-'+ pid}.txt"), val(pid)

  script:
    pid = UUID.randomUUID().toString().substring(0,7)
  
  """
  scanpy-cli pca --n-comps ${n_comps} --svd-solver ${svd_solver} --use-all -f anndata data.h5ad ${ppid + '-' + pid}.h5ad 
  echo "scanpy filter-cells" >> ${ppid +'-'+ pid}.txt
  """
}


// process SCANPY_FILTER_PR {

//   publishDir "$params.output.folder/analysis/scanpy"
//   container "zager/scanpy:3"

//   input:
//     set file("data.h5ad"), val(ppid) from INPUT_CH
//     // each method from params.scanpy.filter_cells
//     // each min_counts from params.scanpy.filter_cells.min_counts
//     // each min_genes from params.scanpy.filter_cells.min_genes
//     // each max_counts from params.scanpy.filter_cells.max_counts
//     // each max_genes from params.scanpy.filter_cells.max_genes

//   output:
//     set file("${ppid + '-' + pid}.h5ad"), file("${ppid +'-'+ pid}.txt"), val(pid) into MONOCLE3_PREPROCESS_CH

//   script:
//     pid = UUID.randomUUID().toString().substring(0,7)

//   """
//   scanpy-cli filter -l -f anndata data.h5ad ${ppid + '-' + pid}.h5ad >> ${ppid +'-'+ pid}.txt
//   // echo "scanpy filter-cells" >> ${ppid +'-'+ pid}.txt
//   """
// }
