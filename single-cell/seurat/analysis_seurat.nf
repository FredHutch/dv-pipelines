#!/usr/bin/env nextflow

INPUT_CH = Channel.from([[
  file(params.output.folder+ '/format/seurat/v3/data.rds'), 
  UUID.randomUUID().toString().substring(0,7)
]])

process SEURAT_FILTER_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), val(ppid) from INPUT_CH
    each low_threshold from params.seurat.filter.low_threshold
    each high_threshold from params.seurat.filter.high_threshold
    each subset_names from params.seurat.filter.subset_names
    each cells_use from params.seurat.filter.cells_use
        
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid) into SEURAT_FILTER_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-filter-cells.R --input-object-file=data.rds --subset-names=${subset_names} --high-thresholds=${high_threshold} --low-thresholds=${low_threshold} --output-object-file=${ppid}-${pid}.rds
  echo "seurat filter-cells subset-names=${subset_names} high-thresholds=${high_threshold} low-thresholds=${low_threshold}" >> ${ppid +'-'+ pid}.txt
  """
}

process SEURAT_NORMALIZE_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid) from SEURAT_FILTER_CH
    each assay_type from params.seurat.normalize.assay_type
    each normalization_method from params.seurat.normalize.normalization_method
    each scale_factor from params.seurat.normalize.scale_factor

  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid) into SEURAT_NORMALIZE_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-normalise-data.R --assay-type=${assay_type} --normalization-method=${normalization_method} --scale-factor=${scale_factor} --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo "seurat normalize-data assay-type=${assay_type} normalization-method=${normalization_method} scale-factor=${scale_factor}" >> ${ppid +'-'+ pid}.txt
  """ 
}

process SEURAT_VARIABLE_GENES_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid) from SEURAT_NORMALIZE_CH
      each dispersion_function from params.seurat.variable_genes.dispersion_function
      each mean_function from params.seurat.variable_genes.mean_function
      each fvg_x_low_cutoff from params.seurat.variable_genes.fvg_x_low_cutoff
      each fvg_y_low_cutoff from params.seurat.variable_genes.fvg_y_low_cutoff
      each fvg_x_high_cutoff from params.seurat.variable_genes.fvg_x_high_cutoff
      each fvg_y_high_cutoff from params.seurat.variable_genes.fvg_y_high_cutoff
    
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid), file("genes.txt") into SEURAT_VARIABLE_GENES_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-find-variable-genes.R --output-text-file=genes.txt --mean-function=${mean_function} --dispersion-function=${dispersion_function} --x-low-cutoff=${fvg_x_low_cutoff} --x-high-cutoff=${fvg_x_high_cutoff} --y-low-cutoff=${fvg_y_low_cutoff} --y-high-cutoff=${fvg_y_high_cutoff} --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo "seurat variable-genes --mean-function=${mean_function} --dispersion-function=${dispersion_function} --x-low-cutoff=${fvg_x_low_cutoff} --x-high-cutoff=${fvg_x_high_cutoff} --y-low-cutoff=${fvg_y_low_cutoff} --y-high-cutoff=${fvg_y_high_cutoff} " >> ${ppid +'-'+ pid}.txt
  """ 
}

process SEURAT_SCALE_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid), file("genes.txt") from SEURAT_VARIABLE_GENES_CH
      each vars_to_regress from params.seurat.scale.vars_to_regress
      each use_umi from params.seurat.scale.use_umi
      each model_use from params.seurat.scale.model_use
      each do_scale from params.seurat.scale.do_scale
      each do_center  from params.seurat.scale.do_center
      each scale_max from params.seurat.scale.scale_max 
      each block_size from params.seurat.scale.block_size 
      each assay_type from params.seurat.scale.assay_type 
      each check_for_norm from params.seurat.scale.check_for_norm
      each min_cells_to_block from params.seurat.scale.min_cells_to_block 
        
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid), file("genes.txt") into SEURAT_SCALE_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-scale-data.R --vars-to-regress=${vars_to_regress} --use-umi=${use_umi} --model-use=${model_use} --do-scale=${do_scale} --do-center=${do_center} --scale-max=${scale_max} --block-size=${block_size} --min-cells-to-block=${min_cells_to_block} --assay-type=${assay_type} --check-for-norm=${check_for_norm} --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo "seurat scale vars-to-regress=${vars_to_regress} use-umi=${use_umi} model-use=${model_use} do-scale=${do_scale} do-center=${do_center} scale-max=${scale_max} block-size=${block_size} min-cells-to-block=${min_cells_to_block} assay-type=${assay_type} check-for-norm=${check_for_norm}" >> ${ppid +'-'+ pid}.txt
  """
}

process SEURAT_PCA_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid), file("genes.txt") from SEURAT_SCALE_CH
      each pcs_compute from params.seurat.pca.pcs_compute
      each use_imputed from params.seurat.pca.use_imputed
    
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid), file("genes.txt") into SEURAT_PCA_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-run-pca.R --pcs-compute=${pcs_compute} --use-imputed=${use_imputed} --output-embeddings-file=embeddings.txt --output-loadings-file=loadings.txt --output-stdev-file=stdev.txt --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo "seurat pca pcs-compute=${pcs_compute} use-imputed=${use_imputed}" >> ${ppid +'-'+ pid}.txt
  """
}

process SEURAT_CLUSTER_PR {

  // publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid), file("genes.txt") from SEURAT_PCA_CH
      each genes_use from params.seurat.cluster.genes_use
      each reduction_type from params.seurat.cluster.reduction_type
      each dims_use from params.seurat.cluster.dims_use
      each k_param from params.seurat.cluster.k_param
      each prune_snn from params.seurat.cluster.prune_snn
      each resolution from params.seurat.cluster.resolution
      each algorithm from params.seurat.cluster.algorithm
    
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid), file("genes.txt") into SEURAT_CLUSTER_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-find-clusters.R --genes-use=${genes_use} --reduction-type=${reduction_type} --dims-use=${dims_use} --k-param=${k_param} --prune-snn=${prune_snn} --resolution=${resolution} --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo "seurat clusters genes-use=${genes_use} reduction-type=${reduction_type} dims-use=${dims_use} k-param=${k_param} prune-snn=${prune_snn} resolution=${resolution}" >> ${ppid +'-'+ pid}.txt
  """
}

process SEURAT_TSNE_PR {

  publishDir "$params.output.folder/analysis/seurat"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("data.rds"), file("ledger.txt"), val(ppid), file("genes.txt") from SEURAT_CLUSTER_CH
    each reduction_use from params.seurat.tsne.reduction_use
    each dims_use from params.seurat.tsne.dims_use
    each do_fast from params.seurat.tsne.do_fast
    
  output:
    set file("${ppid + '-' + pid}.rds"), file("${ppid +'-'+ pid}.txt"), val(pid), file("genes.txt") into SEURAT_TSNE_CH
  
  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  seurat-run-tsne.R --reduction-use=${reduction_use} --dims-use=${dims_use} --do-fast=${do_fast} --input-object-file=data.rds --output-object-file=${ppid}-${pid}.rds
  cat ledger.txt > ${ppid +'-'+ pid}.txt
  echo seurat tsne "seurat reduction-use=${reduction_use} dims-use=${dims_use} do-fast=${do_fast}" >> ${ppid +'-'+ pid}.txt
  """
}