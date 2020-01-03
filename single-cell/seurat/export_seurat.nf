INPUT_CH = Channel
  .fromPath("$params.output.folder/analysis/seurat/*.rds")
  .map { file -> tuple(file.baseName, file) }


process SEURAT_EXTRACT_DIMS_PR {

  publishDir "$params.output.folder/analysis/seurat/$fileName"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set fileName, file("cds.rds") from INPUT_CH

  output:
    set fileName, file("cds.rds"), file("dr_*"), file("attr_col.csv"), file("variance_explained.csv") into SEURAT_EXTRACT_DIMS_CH

  """
  #!/usr/bin/env Rscript
  library(Matrix)
  library(monocle3)

  cds <- readRDS("cds.rds")

  file_prefix <- "${fileName}"

  # reducedDims
  lapply(reducedDimNames(cds), function(x) {
      write.csv(reducedDims(cds)[[x]], row.names=TRUE, quote = FALSE, file = paste0('dr_', tolower(x), ".csv"))
  })

  # clusters/partitions/column attrs
  all_clust <- lapply(names(cds@clusters), function(x) {
      out <- as.data.frame(cbind(cds@clusters[[x]]\$clusters, cds@clusters[[x]]\$partitions))
      names( out) <- paste(names(cds@clusters), c("cluster", "partitions"), sep="_")
      out
  })
  cold <- do.call(cbind, c(as.data.frame(colData(cds)), all_clust))
  write.csv(cold, row.names=TRUE, quote = FALSE, file = paste0("attr_col.csv"))

  # pc variance explained
  pc_var <- data.frame(PC = paste("PC", 1:length(cds@preprocess_aux\$prop_var_expl), sep="_"), Var = cds@preprocess_aux\$prop_var_expl)
  write.csv(pc_var, row.names=FALSE, quote = FALSE, file = paste0("variance_explained.csv"))

  """
}

process SEURAT_EXTRACT_MTXATTR_PR {
  
  publishDir "$params.output.folder/analysis/seurat/$fileName"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set fileName, file("cds.rds"), file("dr_*"), file("attr_col.csv"), file("variance_explained.csv") from SEURAT_EXTRACT_DIMS_CH

  output:
    set fileName, file("cds.rds"), file("attr_col.csv"), file("attr_row.csv"), file("mtx_*") into SEURAT_EXTRACT_MTXATTR_CH

  """
  #!/usr/bin/env Rscript
  library(Matrix)
  library(monocle3)

  cds <- readRDS("cds.rds")

  write.csv(colData(cds), row.names=TRUE, quote = FALSE, file = paste0("attr_col.csv"))
  write.csv(rowData(cds), row.names=TRUE, quote = FALSE, file = paste0("attr_row.csv"))
  count <- t(as(counts(cds), "dgTMatrix"))
  count_df <- data.frame(r = count@i, c = count@j, v = count@x)
  write.csv(count_df, row.names=FALSE, quote = FALSE, file = paste0("mtx_", count@Dim[[1]], '_by_',count@Dim[[2]],'_sparse.csv'))
  """
  
}

process SEURAT_EXTRACT_TRAJECTORY_PR {
  
  errorStrategy 'ignore'  // Would be better to exit detect if umap is in rd before calling

  publishDir "$params.output.folder/analysis/seurat/$fileName"
  container "quay.io/biocontainers/seurat-scripts:0.0.5--r34_1"

  input:
    set fileName, file("cds.rds"), file("attr_col.csv"), file("attr_row.csv"), file("matrix.csv")   from SEURAT_EXTRACT_MTXATTR_CH

  output:
    set fileName, file("cds.rds"), file("trajectory.csv") into SEURAT_EXTRACT_TRAJECTORY_CH

  """
  #!/usr/bin/env Rscript
  library(Matrix)
  library(monocle3)
  library(dplyr)

  cds <- readRDS("cds.rds")
  reduction_method <- "UMAP"
  if(is.null(cds@principal_graph_aux[[reduction_method]]\$dp_mst)) stop("You haven't called learn_graph yet")
  dp_mst <- cds@principal_graph[[reduction_method]]
  if(ncol(reducedDims(cds)[[reduction_method]]) == 2) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]\$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
    edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>%
      dplyr::select_(
      source="sample_name",
      source_prin_graph_dim_1="prin_graph_dim_1",
      source_prin_graph_dim_2="prin_graph_dim_2"),
      by = "source") %>%
    dplyr::left_join(ica_space_df %>%
      dplyr::select_(
      target="sample_name",
      target_prin_graph_dim_1="prin_graph_dim_1",
      target_prin_graph_dim_2="prin_graph_dim_2"),
      by = "target")
  } else {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]\$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2,
      prin_graph_dim_3 = 3) %>%
    dplyr::mutate(sample_name = rownames(.),
      sample_state = rownames(.))
    edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>%
      dplyr::select_(source="sample_name",
      source_prin_graph_dim_1="prin_graph_dim_1",
      source_prin_graph_dim_2="prin_graph_dim_2",
      source_prin_graph_dim_3="prin_graph_dim_3"),
      by = "source") %>%
    dplyr::left_join(ica_space_df %>%
      dplyr::select_(target="sample_name",
      target_prin_graph_dim_1="prin_graph_dim_1",
      target_prin_graph_dim_2="prin_graph_dim_2",
      target_prin_graph_dim_3="prin_graph_dim_3"),
      by = "target")
  }
  write.csv(edge_df, row.names=FALSE, quote = FALSE, file = paste0("trajectory.csv"))
  """
}