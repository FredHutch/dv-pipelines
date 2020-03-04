#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cds_path + '/*_filtered_monocle3_cds.rds')
sample_list = Channel.from(params.input.sample_list)
num_cores = Channel.from(params.input.num_cores)
num_dim = Channel.from(params.monocle.preprocess.num_dim)
norm_method = Channel.from(params.monocle.preprocess.norm_method)
dimred_method = Channel.from(params.monocle.preprocess.reduction_method)
reduce_method = Channel.from(params.monocle.reducedims.reduction_method)
clustering_method = Channel.from(params.monocle.clustercells.clustering_method)
cluster_number = Channel.from(params.monocle.clustercells.cluster_number)

process MON_CLUSTER {
  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/Monocle/Cluster" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'
  label 'gizmo'
  input:
    each cds from cds_path
    val core from num_cores
    each pdim from num_dim
    val pnorm from norm_method
    val preduce from dimred_method
    each rreduce from reduce_method
    each cmethod from clustering_method
    each cnumber from cluster_number

  output:
    file "Monocle_${uuid}.rds" into monocle_cds
    file "Monocle_${uuid}_${rreduce}.png" into monocle_png
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    sample <- metadata(cds)\$Sample
    cds <- preprocess_cds(cds, num_dim = ${pdim}, norm_method = "${pnorm}", method = "${preduce}")
    cds <- reduce_dimension(cds, reduction_method = "${rreduce}", preprocess_method = "${preduce}")
    cds <- cluster_cells(cds, reduction_method = "${rreduce}", k = ${cnumber}, cluster_method = "${cmethod}")
    run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Monocle", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = ${pdim}, preprocess_norm_method = "${pnorm}",
                            preprocess_method = "${preduce}", reduction_method = "${rreduce}", k = ${cnumber}, cluster_method = "${cmethod}", sample = sample)
    metadata(cds) <- list(run_condition = run_condition)
    saveRDS(cds, paste(paste("Monocle", "${uuid}", sep="_"), "rds", sep="."))
    fig <- plot_cells(cds, reduction_method="${rreduce}")
    ggsave(paste(paste("Monocle", "${uuid}", "${rreduce}", sep="_"), "png", sep="."), fig, units = "in", dpi = "retina")
    """

}

process MON_GATHER_METADATA {

  echo false
  publishDir "$params.output.folder/Monocle/Cluster/Metadata" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from monocle_cds.collect()
  output:
    path "Monocle_cluster_metadata.csv" into discard_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  cds_list <- c("${cds_list.join('\",\"')}")
  getClusterSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    cluster_summary <- metadata(cds)\$run_condition
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Monocle_cluster_metadata.csv")
  """
}
