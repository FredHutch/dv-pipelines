#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cds_path)
sample_list = Channel.fromPath(params.input.sample_list)
num_cores = Channel.from(params.input.num_cores)
preprocess_num_dim = Channel.from(params.preprocess.num_dim)
preprocess_norm_method = Channel.from(params.preprocess.norm_method)
preprocess_reduce_method = Channel.from(params.preprocess.reduction_method)

INPUT_CH = Channel.from([[
  path(params.input.cds_path), 
  val(params.input.sample_list), 
  val(params.input.row), 
  UUID.randomUUID().toString().substring(0,7)
]])



process MON_CLUSTER {
  echo false
  publishDir "$params.output.folder/Monocle/Preprocess/CDS"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    path cds from cds_path
    each sample from sample_list
    val core from num_cores
    each pre_dim from preprocess_num_dim
    val pre_norm from preprocess_norm_method
    val pre_reduce from preprocess_reduction_method
    each red_reduce from reducedims_reduce_method

  output:
    path "${study}_GEX_h5_samplesheet.csv" into gex_h5sheet
    path "${study}_GEX_samplesheet.csv" into gex_samplesheet
    path "${study}_VDJ_samplesheet.csv" into vdj_samplesheet
    path "${study}_VDJ_analysis_samplesheet.csv" into vdj_analysissheet
        
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    gene_metadata <- rowData(cds)
    cell_metadata <- colData(cds)
    row.names(cell_metadata) <- cell_metadata\$Barcode
    cds <- new_cell_data_set(expression_data = counts(cds), gene_metadata = gene_metadata, cell_metadata = cell_metadata)
    cds <- preprocess_cds(cds, num_dim=${pre_dim}, norm_method =${pre_norm}, method = ${pre_reduce})
    cds <- reduce_dimensions(cds, reduction_method=${red_reduce}, preprocess_method=${pre_reduce})

    
    """

}

