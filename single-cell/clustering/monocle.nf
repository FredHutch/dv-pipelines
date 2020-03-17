#!/usr/bin/env nextflow
// Create input channels
cds_path = Channel.fromPath(params.input.cds_path + '/*_monocle3_cds.rds')
study_name = Channel.from(params.input.study_name)
num_dim = Channel.from(params.monocle.preprocess.num_dim)
norm_method = Channel.from(params.monocle.preprocess.norm_method)
dimred_method = Channel.from(params.monocle.preprocess.reduction_method)
reduce_method = Channel.from(params.monocle.reducedims.reduction_method)
clustering_method = Channel.from(params.monocle.clustercells.clustering_method)
cluster_number = Channel.from(params.monocle.clustercells.cluster_number)
align_method = Channel.from(params.monocle.aligncds.preprocess_method)
align_dims = Channel.from(params.monocle.aligncds.num_dim)
align_group = Channel.from(params.monocle.aligncds.alignment_group)
//Duplicate input channels 
cds_path.into{mon_combine; mon_data}
//Combine CDS object and run monocle cluster
process MON_COMBINE {
  echo false
  scratch "$task.scratch"
  module "R/3.6.1-foss-2016b-fh2"
  label "gizmo_largenode"
  
  input:
    val cds_list from mon_combine.collect()
    val study_name from study_name
  output:
    path "${study_name}_aggregated.rds" into mon_aggr

  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    library(purrr)
    set.seed(12357)

    cds_files <- c("${cds_list.join('\",\"')}")
    cds_list <- cds_files %>% map(readRDS)
    merged <- combine_cds(cds_list)
    metadata(merged)\$Sample <- "${study_name}_aggregated"

    #Combine metdata
    combine_metrics <- function(cds) {
      sample <- metadata(cds)\$Sample
      metrics <- metadata(cds)\$perCellQCMetrics_filtered %>% as_tibble() %>% mutate(Sample = sample)
      return(metrics)
    }
    combine_info <- function(cds) {
      sample <- metadata(cds)\$Sample
      info <- metadata(cds)\$study_info
      return(info)
    }
    combine_vdj <- function(cds) {
      sample <- metadata(cds)\$Sample
      if (metadata(cds)\$vdj_raw != NULL) {
        metadata(cds)\$vdj_raw <- metadata(cds)\$vdj_raw %>% mutate(Sample = sample)
      }
      info <- metadata(cds)\$study_info
      return(info)
    }
    metadata(merged)\$perCellQCMetrics <- cds_list %>% map(combine_metrics) %>% purrr::reduce(rbind)
    metadata(merged)\$study_info <- cds_list %>% map(combine_info) %>% purrr::reduce(rbind)
    metadata(merged)\$vdj_raw <- cds_list %>% map(combine_vdj) %>% purrr::reduce(rbind)        
    saveRDS(merged, "${study_name}_aggregated.rds")
    """
}
//Combine aggregated and sample channels
mon_path = mon_data.mix( mon_aggr )
//Preprocess cells
process MON_PREPROCESS {
  echo false
  scratch "$task.scratch"
  module "R/3.6.1-foss-2016b-fh2"
  label "gizmo_largenode"
  
  input:
    each cds from mon_path
    each pdim from num_dim 
    each pnorm from norm_method
    each preduce from dimred_method

  output:
    path filename into mon_prep

  script:
    filename = cds.getName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    sample <- metadata(cds)\$Sample
    cds <- preprocess_cds(cds, num_dim=${pdim}, norm_method="${pnorm}", method="${preduce}")
    metadata(cds)\$clustering_params <- tibble(initial_reduction = c("${preduce}"), 
                                               intial_dimension = c(${pdim}),
                                               normalization_method = c("${pnorm}"))
    saveRDS(cds, "${filename}")
    """
}
// Select merged sample from pre-process channel
mon_merged = Channel.create()
mon_samples = Channel.create()
mon_prep.choice(mon_merged, mon_samples) { a -> a =~ /^.*_aggregated.rds/ ? 0 : 1}
// Batch correct merged CDS object
process MON_BATCHCORRECT {
  echo false
  scratch "$task.scratch"
  module "R/3.6.1-foss-2016b-fh2"
  label "gizmo_largenode"
  
  input:
    each cds from mon_merged
    val pmethod from align_method
    val ndim from align_dims 
    val agroup from align_group

  output:
    path "${filename}" into mon_aligned

  script:
    filename = cds.getName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cds <- align_cds(cds = cds, num_dim = ${ndim}, 
                     alignment_group = "${agroup}", preprocess_method = "${pmethod}")
    saveRDS(cds, "${filename}")
    """
}
//Combine channels again
mon_recombine = mon_samples.mix(mon_aligned)
// Reduce dimensions
process MON_REDUCEDIMS {
  echo false
  scratch "$task.scratch"
  module "R/3.6.1-foss-2016b-fh2"
  label "gizmo_largenode"
  
  input:
    each cds from mon_recombine
    each rreduce from reduce_method

  output:
    path filename into mon_reduce 

  script:
    filename = cds.getName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    preprocess_method <- metadata(cds)\$clustering_params\$initial_reduction[1]
    cds <- reduce_dimension(cds, reduction_method = "${rreduce}", 
                            preprocess_method = preprocess_method, cores=$task.cpus)
    metadata(cds)\$clustering_params <- metadata(cds)\$clustering_params %>% add_column(reduction_method = c("${rreduce}"))
    saveRDS(cds, "${filename}")
    """
}
//CLuster cells
process MON_CLUSTER {
  echo false
  scratch "$task.scratch"
  module "R/3.6.1-foss-2016b-fh2"
  label "gizmo_largenode"
  
  input:
    each cds from mon_reduce
    each cmethod from clustering_method
    each cnumber from cluster_number

  output:
    path "${filename}_${uuid}.rds" into mon_cluster

  script:
    filename = cds.getSimpleName() - ~/_filtered_monocle3_cds$/
    uuid = UUID.randomUUID().toString().substring(0,7)   
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    sample <- metadata(cds)\$Sample
    reduction_method <- metadata(cds)\$clustering_params\$reduction_method[1]
    cds <- cluster_cells(cds, reduction_method = reduction_method, k = ${cnumber}, cluster_method = "${cmethod}")
    metadata(cds)\$clustering_params <- metadata(cds)\$clustering_params %>% add_column(nearest_neighbors = c(${cnumber}), 
                                        clustering_method = c("${cmethod}"), uuid = c("${uuid}"))
    saveRDS(cds, "${filename}_${uuid}.rds")
    """
}
//Duplicate CDS channel
mon_cluster.into{mon_plot; mon_gather; mon_write}
//Cobine CDS objects and normalization
//Plot cluster scatter plot in UMAP/tSNE dimensions
process MON_PLOT {
  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/Monocle/Cluster/Figures" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'
  label 'gizmo'
  input:
    each cds from mon_plot

  output:
    file "${sample_name}_*.png" into monocle_png
  script:
    sample_name = cds.getSimpleName() - ~/_\w{7}$/
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(monocle3)
    set.seed(12357)

    cds <- readRDS("${cds}")
    reduction_method <- metadata(cds)\$clustering_params\$reduction_method[1]
    uuid <- metadata(cds)\$clustering_params\$uuid[1]
    sample = metadata(cds)\$Sample 
    if (str_detect(sample, "aggregated")) {
      fig <- plot_cells(cds, reduction_method=reduction_method, color_cells_by="Sample")  
    } else {
      fig <- plot_cells(cds, reduction_method=reduction_method)
    }
    ggsave(paste(paste("${sample_name}", uuid, sep="_"), "png", sep="."), fig, 
                 units = "in", height = 10, width = 15, dpi = "retina")
    """

}
// Gather metadata information from all samples
process MON_GATHER {
  echo false
  publishDir "$params.output.folder/Monocle/Cluster/Metadata" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from mon_gather.collect()
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
    cluster_summary <- metadata(cds)\$clustering_params
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Monocle_cluster_metadata.csv")
  """
}
// Write CDS object
process MON_WRITE {
  echo false
  publishDir "$params.output.folder/Monocle/Cluster/CDS" , mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'
  label 'gizmo'

  input:
    each cds from mon_write
  output:
    path "${filename}" into cds_output
  
 script:
    filename = cds.getName()
    """
    #!/usr/bin/env Rscript
    library(monocle3)
    cds <- readRDS("${cds}")
    uuid <- metadata(cds)\$clustering_params\$uuid[1]
    saveRDS(cds, "${filename}")
    """
}
