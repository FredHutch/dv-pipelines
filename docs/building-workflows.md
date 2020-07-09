# So You Want to Build a Workflow

This is documentation on how to develop a workflow, or convert an existing shell script to run as a workflow.

There are a few main things to consider:

## Important Variables

* Storage (S3 location)
* 

* [Input records + ](inputs.md)
* Parameters JSON blob
* Batch.json
* 






The Fred Hutch Data Visualization Center is developing a platform to run and visualize genomic analysis online.
The platform is composed of multiple projects including:

1. A Reference Data API - Containing 2TB of information on genes and variants from over 25 public repositories.
2. A Web Visualization Components - That contain genomic specific visualization built with WebGL to facilitate the rendering of multi-million point datasets.
3. A Web Framework / CMS - To facilitate the aggregation of visualization components into "Atlas" websites.
4. Cloud Infrastructure to support the storage and retrieval of data through rest, graphql and websockets.
5. Data Visualization Pipelines - Which allows for community submission of compatible pipelines.

This repository contains a collection of single cell pipelines compatible with the Fred Hutch Visualization Center Platform

## Repo Contents

This repository contains a Platform compatible pipelines. Our inital focus will be to create scRNA-Seq and ImmunoSeq analysis.

The repository contains nextflow scripts to run different workflows for scRNA-Seq analysis as well as MiXCR workflow for AIRR-Seq data analysis. The scRNA-Seq workflows are divide into five logical steps.  Over time we will be working to convert all of these scripts to be compatible with our platform.  The "example" files provide a basic template for contributing pipelines.

1. Count matrix generation
2. Preprocessing
3. Clustering
4. Cell-type classification
5. Differential expression analysis
6. Trajectory analysis






## Batch Call To Trigger Example

aws batch submit-job \
     --job-name example-job \
     --job-queue general-purpose-queue \
     --job-definition nextflow \
     --container-overrides command=FredHutch/dv-pipelines/example.nf, \
        "--dataset", "f49a1d9f-a165-4b0d-a53a-7adab578d799"

[Inputs (forms + records)](inputs.md)
[Workflows](workflows.md)
[Visualizations](visualizations.md)


