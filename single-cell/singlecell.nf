#!/usr/bin/env nextflow

INPUT_CH = Channel.from([[
  file(params.input.exp), 
  file(params.input.col), 
  file(params.input.row), 
  UUID.randomUUID().toString().substring(0,7)
]])


// Reacts To Input Channel - Outputs 
INPUT_CH.into {
  TENX_INPUT_CH
  CDS_INPUT_CH
}

process TENX_LOAD_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file('exp.rds'), file('col.rds'), file('row.rds'), val(ppid) from TENX_INPUT_CH

  output:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv'), val(pid) into TENX_LOAD_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

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

process SCANPY_LOAD_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/scanpy-scripts:0.0.3--py37_1"

  input:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv'), val(ppid) from SCANPY_INPUT_CH

  output:
    set file("${ppid + '-' + pid}.h5ad"), file('ledger.txt'), val(pid)

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
    scanpy-read-10x.py -d \$PWD/ -o ${ppid + '-' + pid}.h5ad -F anndata
    echo "$ppid-$pid scanpy create" >> ledger.txt
  """
}


process SEURAT_LOAD_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file('matrix.mtx'), file('barcodes.tsv'), file('cells.tsv'), file('genes.tsv'), val(ppid) from SEURAT_INPUT_CH

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_LOAD_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  Rscript /usr/local/bin/seurat-read-10x.R -d \$PWD -o seurat-tmp.rds
  Rscript /usr/local/bin/seurat-create-seurat-object.R -i seurat-tmp.rds -o ${ppid + '-' + pid}.rds
  echo "${ppid + '-' + pid} seurat create" >> ledger.txt
  """
}


process SEURAT_FILTERCELLS_PR {

  errorStrategy 'ignore'  // Would be better to conditionally broacast if MONOCLE_PARTITION_CH Reduction Method = UMAP

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_LOAD_CH
    each low_threshold from params.seurat.filter_cells.low_threshold
    each high_threshold from params.seurat.filter_cells.high_threshold

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_FILTERCELLS_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

  """
  Rscript /usr/local/bin/seurat-filter-cells.R -i seurat.rds -o ${ppid + '-' + pid}.rds --low-thresholds=${low_threshold}  --high-thresholds=${high_threshold}
  echo "${ppid + '-' + pid} seurat filtercells low_threshold=${low_threshold} high_threshold=${high_threshold}" >> ledger.txt
  """
}

process SEURAT_NORMALIZE_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"
  
  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_FILTERCELLS_CH
    each assay_type from params.seurat.normalize.assay_type
    each normalization_method from params.seurat.normalize.normalization_method
    each scale_factor from params.seurat.normalize.scale_factor

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_NORMALIZE_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

  """
  Rscript /usr/local/bin/seurat-normalise-data.R -i seurat.rds -o ${ppid + '-' + pid}.rds -a ${assay_type} -n ${normalization_method} -s ${scale_factor} 
  echo "${ppid + '-' + pid} seurat normalize assay-type=${assay_type} normalization_method=${normalization_method} scale_factor=${scale_factor}" >> ledger.txt
  """
}

process SEURAT_VARIABLEGENES_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"
  
  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_NORMALIZE_CH
    each mean_function from params.seurat.variable_genes.mean_function
    each dispersion_function from params.seurat.variable_genes.dispersion_function
    each fvg_x_low_cutoff from params.seurat.variable_genes.fvg_x_low_cutoff
    each fvg_y_low_cutoff from params.seurat.variable_genes.fvg_y_low_cutoff
    each fvg_x_high_cutoff from params.seurat.variable_genes.fvg_x_high_cutoff
    each fvg_y_high_cutoff from params.seurat.variable_genes.fvg_y_high_cutoff

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_VARIABLEGENES_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

  """
  Rscript /usr/local/bin/seurat-find-variable-genes.R -i seurat.rds -m ${mean_function} -d ${dispersion_function} -l ${fvg_x_low_cutoff} -j ${fvg_x_high_cutoff} -y ${fvg_y_low_cutoff} -z ${fvg_y_high_cutoff} -o ${ppid + '-' + pid}.rds -t ${ppid + '-' + pid}.txt
  echo "$ppid-$pid seurat variableGenes mean-function=${mean_function} dispersion-function=${dispersion_function} fvg-x-low-cutoff=${fvg_x_low_cutoff} fvg-x-high-cutoff=${fvg_x_high_cutoff} fvg-y-low-cutoff=${fvg_y_low_cutoff} fvg-y-high-cutoff=${fvg_y_high_cutoff}" >> ledger.txt
  """
}

process SEURAT_SCALEDATA_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"
  
  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_VARIABLEGENES_CH
    each model_use  from params.seurat.scale_data.model_use
    each do_scale from params.seurat.scale_data.do_scale
    each do_center from params.seurat.scale_data.do_center 
    each scale_max from params.seurat.scale_data.scale_max
    each block_size from params.seurat.scale_data.block_size
    each min_cells_to_block from params.seurat.scale_data.min_cells_to_block
    each assay_type from params.seurat.scale_data.assay_type
    each check_for_norm from params.seurat.scale_data.check_for_norm

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_SCALEDATA_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

  """
  Rscript /usr/local/bin/seurat-scale-data.R -i seurat.rds -o ${ppid + '-' + pid}.rds -m ${model_use} -s ${do_scale} -c ${do_center} -x ${scale_max} -b ${block_size} -d ${min_cells_to_block} -a ${assay_type} -n ${check_for_norm}
  echo "$ppid-$pid seurat scaleData m ${model_use} do_scale=${do_scale} do_center=${do_center} scale_max=${scale_max} block_size=${block_size} min_cells_to_block=${min_cells_to_block} assay_type=${assay_type} check_for_norm=${check_for_norm}" >> ledger.txt
  """
}

process SEURAT_PCA_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"
  
  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_SCALEDATA_CH
    each pcs_compute from params.seurat.pca.pcs_compute
    each use_imputed from params.seurat.pca.use_imputed

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_PCA_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

    """
    Rscript /usr/local/bin/seurat-run-pca.R -i seurat.rds -o ${ppid + '-' + pid}.rds -p ${pcs_compute} -m ${use_imputed}
    echo "$ppid-$pid seurat pca pcs_compute=${pcs_compute} use_imputed=${use_imputed}" >> ledger.txt
    """
}

process SEURAT_TSNE_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"
  
  input:
    set file("seurat.rds"), file("ledger.txt"), val(ppid) from SEURAT_PCA_CH
    each reduction_use from params.seurat.tsne.reduction_use
    each do_fast from params.seurat.tsne.do_fast

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into SEURAT_TSNE_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

    """
    Rscript /usr/local/bin/seurat-run-tsne.R -i seurat.rds -o ${ppid + '-' + pid}.rds -f ${do_fast} -r ${reduction_use}
    echo "$ppid-$pid seurat tsne  do_fast=${do_fast} reduction_use=${reduction_use}" >> ledger.txt
    """
}

process CDS_LOAD_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file('exp.rds'), file('col.rds'), file('row.rds'), val(ppid) from CDS_INPUT_CH

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into CDS_LOAD_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7) 

  """
  monocle3 create -F cds3 --cell-metadata=col.rds --gene-annotation=row.rds --expression-matrix=exp.rds ${ppid + '-' + pid}.rds
  echo "$ppid-$pid monocle3 create" >> ledger.txt
  """
}

process MONOCLE3_PREPROCESS_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from CDS_LOAD_CH
    each method from params.monocle3.preprocess.method
    each num_dim from params.monocle3.preprocess.num_dim
    each norm_method from params.monocle3.preprocess.norm_method
    each pseudo_count from params.monocle3.preprocess.pseudo_count

  output:
    set file("${ppid + '-' + pid}.rds"), file('ledger.txt'), val(pid) into MONOCLE3_PREPROCESS_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  monocle3 preprocess -f cds3 -F cds3 --method=${method} --num-dim=${num_dim} --norm-method=${norm_method} --pseudo-count=${pseudo_count} cds.rds ${ppid}-${pid}.rds
  echo "$ppid-$pid monocle3 preprocess method=${method} num-dim=${num_dim} norm-method=${norm_method} pseudo-count=${pseudo_count}" >> ledger.txt
  """
}

process MONOCLE3_REDUCEDIM_PR {

  publishDir "$params.output.folder"
  container 'quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1'

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from MONOCLE3_PREPROCESS_CH
    each max_components from params.monocle3.reduce_dim.max_components
    each steps from params.monocle3.reduce_dim.steps

  output:
    set file("${ppid + '-'+ pid}.rds"), file("ledger.txt"), val(pid) into MONOCLE3_REDUCEDIM_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  preprocess=\$(echo ${steps} | cut -f1 -d-)
  reduction=\$(echo ${steps} | cut -f2 -d-)
  echo "$ppid-$pid monocle3 reduceDim preprocess-method=\$preprocess reduction-method=\$reduction max-components=${max_components}" >> ledger.txt
  monocle3 reduceDim -f cds3 -F cds3 --preprocess-method=\$preprocess --reduction-method=\$reduction --max-components=${max_components} cds.rds ${ppid}-${pid}.rds
  """
}

process MONOCLE3_PARTITION_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from MONOCLE3_REDUCEDIM_CH
    each reduction_method from params.monocle3.partition.reduction_method
    each knn from params.monocle3.partition.knn
    each louvain_iter from params.monocle3.partition.louvain_iter
    each partition_qval from params.monocle3.partition.partition_qval
  output:
    set file("${ppid + '-'+ pid}.rds"), file("ledger.txt"), val(pid) into MONOCLE3_PARTITION_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  monocle3 partition -f cds3 -F cds3 --reduction-method=${reduction_method} --knn=${knn} --louvain-iter=${louvain_iter} cds.rds ${ppid+'-'+pid}.rds
  echo "$ppid-$pid monocle3 partition reduction-method=${reduction_method} knn=${knn} louvain-iter=${louvain_iter}" >> ledger.txt
  """
}

process MONOCLE3_LEARNGRAPH_PR {

  errorStrategy 'ignore'  // Would be better to conditionally broacast if MONOCLE_PARTITION_CH Reduction Method = UMAP

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from MONOCLE3_PARTITION_CH
    each euclidean_distance_ratio from params.monocle3.learn_graph.euclidean_distance_ratio
    each geodesic_distance_ratio from params.monocle3.learn_graph.geodesic_distance_ratio
    each minimal_branch_len from params.monocle3.learn_graph.minimal_branch_len

  output:
    set file("${ppid + '-'+ pid}.rds"), file("ledger.txt"), val(pid) into MONOCLE3_LEARNGRAPH_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  monocle3 learnGraph -f cds3 -F cds3 --euclidean-distance-ratio=${euclidean_distance_ratio} --geodesic-distance-ratio=${geodesic_distance_ratio} --minimal-branch-len=${minimal_branch_len} cds.rds ${ppid}-${pid}.rds
  echo "$ppid-$pid monocle3 learnGraph euclidean-distance-ratio=${euclidean_distance_ratio} geodesic-distance-ratio=${geodesic_distance_ratio} minimal-branch-len=${minimal_branch_len}" >> ledger.txt
  """
}


process MONOCLE3_ORDERCELLS_PR {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from MONOCLE3_LEARNGRAPH_CH
    each cell_phenotype from params.monocle3.order_cells.cell_phenotype
    each reduction_method from params.monocle3.order_cells.reduction_method
  
  output:
    set file("${ppid + '-'+ pid}.rds"), file("ledger.txt"), val(pid) into MONOCLE_ORDERCELLS_CH

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  echo "$ppid-$pid monocle3 orderCells cell-phenotype=${cell_phenotype} reduction-method=${reduction_method}" >> ledger.txt
  monocle3 orderCells -f cds3 -F cds3 --cell-phenotype=${cell_phenotype} --reduction-method=${reduction_method} --root-type=NULL cds.rds ${ppid}-${pid}.rds
  """
}
/*
process monocle3_diffexp {

  publishDir "$params.output.folder"
  container "quay.io/biocontainers/monocle3-cli:0.0.3--py37r36hc9558a2_1"

  input:
    set file("cds.rds"), file("ledger.txt"), val(ppid) from MONOCLE_ORDERCELLS_CH
    each neighbor_graph from params.monocle3.diff_exp.neighbor_graph
    each knn from params.monocle3.diff_exp.knn
    each alternative from params.monocle3.diff_exp.alternative

  output:
    set file("${ppid + '-'+ pid}.rds"), file("ledger.txt"), val(pid)

  script:
    pid = UUID.randomUUID().toString().substring(0,7)

  """
  echo "$ppid-$pid monocle3 diffExp neighbor-graph=${neighbor_graph} reduction-method=${params.monocle3.order_cells.reduction_method} knn=${knn} method=${params.monocle3.diffexp.method} alternative=${alternative}" >> ledger.txt
  monocle3 diffExp -f cds3 -F cds3 --neighbor-graph=${neighbor_graph} --reduction-method=${params.monocle3.order_cells.reduction_method} --knn=${knn} --method=${ params.monocle3.diffexp.method} --alternative=${alternative} cds.rds ${ppid}-${pid}.rds
  """
}
*/
