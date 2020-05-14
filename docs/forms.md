# Forms

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


### scRNA-Seq and AIRR-Seq tools incorporated

1. MiXCR
2. Monocle3
3. Seurat
4. Scanpy
5. Scater
6. Scran
7. DropletUtils

## Batch Call To Trigger Example

aws batch submit-job \
     --job-name example-job \
     --job-queue general-purpose-queue \
     --job-definition nextflow \
     --container-overrides command=FredHutch/dv-pipelines/example.nf, \
        "--dataset", "f49a1d9f-a165-4b0d-a53a-7adab578d799"

## Example Files Required To Register Analysis

Below is a list of files that must exist in the root of this repo for to register the analysis with the system.
Files should be prefixed with the technology.phase.library.version.  For example "singlecell.clustering.monocle.3"

For detailed information on how to format json object refer to:

Documentation Site: <https://react-jsonschema-form.readthedocs.io/en/latest/advanced-customization/>

Online Test Harness: <https://cybertec-postgresql.github.io/rjsf-material-ui/>

### File: example.form.json

This file contains a json object that is used to describe the inputs to the nextflow script.  
It's primary purpose is to render a webform to collect parameters from an end user.  
A secondary benifit is it can be used to validate input.

```JSON
{
  "title": "Test Script",
  "description": "Tests Basic ScanPY Function",
  "type": "object",
  "required": [
    "n_comps",
    "zero_center"
  ],
  "properties": {
    "n_comps": {
      "type": "integer",
      "title": "Number of Components"
    },
    "zero_center": {
      "type": "boolean",
      "title": "Zero Center"
    }
  }
}
```

### File: example.form.ui.json

This file can used to further customize the form.  It can include help text, enumeration values, choice of UI Components
If you are comfortable with a basic form then this file can simply include an empty JSON object (EG {})

```JSON
{
  "n_comps": {
    "ui:widget": "updown"
  }
}
```

### File: example.form.values.json

This file contains a json object with the default values to render in the form.
When the user specifies parameters to run they will be serialized in the same format + placed in the nextflow working directory as
"params.json"

```JSON
{
  "zero_center": false,
  "n_comps": 50
}
```

