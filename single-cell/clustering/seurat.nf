#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cds_path + "/*_seurat_cds.rds")
sample_list = Channel.fromPath(params.input.sample_list)
nmethods = Channel.from(params.seurat.normalize.method)
smethods = Channel.from(params.seurat.variable.method)
nfeatures = Channel.from(params.seurat.variable.nfeatures)
ndims = Channel.from(params.seurat.neighbors.dims)
ncenters = Channel.from(params.seurat.neighbors.centers)
algorithms = Channel.from(params.seurat.cluster.algorithm)
resolutions = Channel.from(params.seurat.cluster.resolution)

process SEURAT_CLUSTER {
  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/Seurat/Cluster/CDS", mode: 'copy'
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.6.7-foss-2016b-fh2'
  label 'gizmo'
  input:
    each cds from cds_path
    each nmethod from nmethods
    each smethod from smethods
    each nfeature from nfeatures
    each ndim from ndims
    each ncenter from ncenters
    each algorithm from algorithms
    each resolution from resolutions
    

  output:
    file "Seurat_${uuid}.rds" into seurat_cds
    file "Seurat_${uuid}_umap.png" into seurat_umaps
    file "Seurat_${uuid}_tsne.png" into seurat_tsnes
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)    
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(scater)
    library(Seurat)
    library(tidyverse)
    set.seed(12357)

    cds <- readRDS("${cds}")
    #Normalize the data
    cds <- NormalizeData(cds, normalization.method = "${nmethod}", scale.factor = 10000)
    #Find genes that are show highly variable expression and consider a subset of these for further analysis
    cds <- FindVariableFeatures(cds, selection.method = "${smethod}", nfeatures = ${nfeature})
    gene_names <- rownames(cds)
    # Scale the data so as to bring the mean to 0 and variance to 1, this prevents the greater influence of high expression values
    cds <- ScaleData(cds, features = gene_names)
    # Run PCA
    cds <- RunPCA(cds, features = VariableFeatures( object = cds), npcs=${ndim})
    #Cluster cells
    cds <- FindNeighbors(cds, dim = 1:${ndim}, k=${ncenter})
    cds <- FindClusters(cds, resolution = ${resolution}, algorithm = ${algorithm})
    cds <- RunUMAP(cds, dims = 1:${ndim})
    cds <- RunTSNE(cds, dims = 1:${ndim})
    run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Seurat", "${uuid}", sep="_"), "rds", sep="."), 
			   normalize_nmethod = "${nmethod}", variable_exp_method = "${smethod}",
                            variable_exp_features = ${nfeature}, pca_dimensions = ${ndim}, 
			k = ${ncenter}, resolution = ${resolution}, cluster_method = ${algorithm}, sample = sample)
    cds_sce <- as.SingleCellExperiment(cds)
    metadata(cds_sce) <- list(seurat_run_condition = run_condition)
    saveRDS(cds_sce, paste(paste("Seurat", "${uuid}", sep="_"), "rds", sep="."))
    umap_fig <- DimPlot(cds, reduction = "umap")
    ggsave(paste(paste("Seurat", "${uuid}", "umap", sep="_"), "png", sep="."), umap_fig, units = "in", dpi = "retina")
    tsne_fig <- DimPlot(cds, reduction = "tsne")
    ggsave(paste(paste("Seurat", "${uuid}", "tsne", sep="_"), "png", sep="."), tsne_fig, units = "in", dpi = "retina")
    """

}

process Seurat_GATHER_METADATA {

  echo false
  publishDir "$params.output.folder/Seurat/Cluster" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from seurat_cds.collect()
  output:
    path "Seurat_cluster_metadata.csv" into discard_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(Seurat)
  library(DropletUtils)
  set.seed(12357)
  cds_list <- c("${cds_list.join('\",\"')}")
  getClusterSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    cluster_summary <- metadata(cds)\$seurat_run_condition
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Seurat_cluster_metadata.csv")
  """
}
