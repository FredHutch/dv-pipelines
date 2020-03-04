#!/usr/bin/env nextflow

cds_path = Channel.from(params.input.cds_path)
sample_list = Channel.from(params.input.sample_list)
nfeatures = Channel.from(params.scanpy.variable.nfeatures)
ndims = Channel.from(params.scanpy.pca.dims)
ncenters = Channel.from(params.scanpy.neighbors.centers)
nmethods = Channel.from(params.scanpy.neighbors.method)
resolutions = Channel.from(params.scanpy.cluster.resolution)

process SCANPY_CLUSTER {
  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/Scanpy/Cluster/CDS", mode : 'copy'
  module 'Python/3.7.4-foss-2016b'
  input:
    val cds from cds_path
    each sample from sample_list
    each nfeature from nfeatures
    each ndim from ndims
    each ncenter from ncenters
    each nmethod from nmethods
    each resolution from resolutions

  output:
    file "Scanpy_${uuid}.h5ad" into scanpy_cds
    file "Scanpy_${uuid}_UMAP.png" into scanpy_umaps
    file "Scanpy_${uuid}_tSNE.png" into scanpy_tsne
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)            
  
    """
    #!/usr/bin/env python3
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from matplotlib import pyplot as plt
    np.random.seed(12357)

    cds = sc.read_10x_mtx("${cds}/${sample}_filtered_matrix", var_names = "gene_symbols")
    uuid = "${uuid}" #"${UUID.randomUUID().toString().substring(0,7)}"
    mito_genes = cds.var_names.str.startswith("MT-")
    cds.obs['percent_mito'] = np.sum(cds[:, mito_genes].X, axis =1).A1/ np.sum(cds.X, axis=1).A1
    cds.obs['n_counts'] = cds.X.sum(axis=1).A1
    #Normalize the data
    cds = sc.pp.normalze_total(cds, target_sum=1e4)
    cds = sc.pp.log1p(cds)
    cds.raw = cds
    sc.pp.highly_variable_genes(cds, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes= ${nfeature})
    cds = cds[:, cds.var.highly_variable]
    sc.pp.regress_out(cds, ['n_counts', 'percent_mito'])
    sc.pp.scale(cds, max_value=10)
    # PCA analysis
    sc.tl.pca(cds, n_comps=${ndim})
    # Generate neighbor graph 
    sc.pp.neighbors(cds, n_neighbors=${ncenter}, n_pcs=${ndims}, method='${nmethod}')
    # Generate UMAP and tSNE
    sc.tl.tsne(cds, n_pcs=${ndim})
    sc.tl.umap(cds)
    sc.tl.leiden(cds, resolution=${resolution})
    sc.tl.louvain(cds, resolution=${resolution})
    out_file = 'Scanpy_{0}.h5ad'.format("${uuid}")
    run_conditions = {'uuid': "${uuid}", 'outfile': out_file, 'sample': '${sample}', 'variable_gene_features': ${nfeature}, 'pca_dimensions': ${ndim}, 'k': ${ncenter}, 'neighbor_joining_method': '${nmethod}', 'resolution': ${resolution}}
    cds.uns['run_conditions'] = run_conditions
    cds.write(out_file)
    tsne_plot = sc.pl.tsne(cds)
    umap_plot = sc.pl.umap(cds)
    umap_plot.savefig("Scanpy_{0}_UMAP.png".format("${uuid}"))
    tsne_plot.savefig("Scanpy_{0}_tSNE.png".format("${uuid}"))
    """
}
