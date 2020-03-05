# Overview

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

### Test Dynamo DB Records

#### Dataset (Test)

```JSON
{
  "id": "f49a1d9f-a165-4b0d-a53a-7adab578d799",
  "configuration": "bebd6d26-59d8-481e-b6b7-99a5823acbfa",
  "datasets": ["8ca1f73e-b7e5-4b66-84b4-ca11ad14c528"],
  "name": "Loom Test Analysis",
  "desc": "Loom Test Analysis of Census of Immune Cells",
  "groupsEdit": ["Admin"],
  "groupsView": ["Admin"],
  "s3path": "S3://dvc-portal-data/hca/f49a1d9f-a165-4b0d-a53a-7adab578d799",
  "status": "success"
}
```

#### Configuration (Test)

```JSON
{
  "id": "bebd6d26-59d8-481e-b6b7-99a5823acbfa",
  "analysis": "9ccc42ac-e566-4f07-9ca2-135ac4ec2bb0",
  "desc": "Human Cell Atlas Inport",
  "groupsEdit": ["Admin"],
  "groupsView": ["Admin"],
  "info": {
    "date": 20200303,
    "source": "data.humancellatlas.org"
  },
  "name": "Human Cell Atlas Inport",
  "parameters": {
      "zero_center": false,
      "n_comps": 50
  },
  "type": "one-time"
}
```

#### Analysis (Test)

```JSON
{
  "id": "9ccc42ac-e566-4f07-9ca2-135ac4ec2bb0",
  "analysis": [],
  "name": "Test Analysis",
  "desc": "Basic Loom Test Analysis",
  "info": {},
  "runtime": "NEXTFLOW",
  "executable": "https://github.com/FredHutch/dv-pipelines/example.nf",
  "schema": {
    "form": {
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
    },
    "ui": {
      "n_comps": {
        "ui:widget": "updown"
      }
    },
    "values": {
      "zero_center": false,
      "n_comps": 50
    }
  },
  "files": ["data.loom"],
  "groupsView": ["Admin"],
  "groupsEdit": ["Admin"]
}
```

### Seed DynamoDB Records

#### Dataset (Seed)

```JSON
{
  "id": "8ca1f73e-b7e5-4b66-84b4-ca11ad14c528",
  "configuration": "fa16a5e9-9c96-4ffb-b5ac-26eecc4b6eb7",
  "datasets": [],
  "desc": "Diverse cells of the immune system maintain and protect tissue function, integrity, and homeostasis upon changes in functional demands and diverse perturbations. Recent advances such as massively parallel single-cell RNA-sequencing and sophisticated computational methods help shed new light on this complexity. This immune cell census aims to profile up to 2M immunocytes, the first tranche of this is currently available. With computational methods optimized to a massive scale, we can readily identify cell types and markers, as well as the process of hematopoietic differentiation. The high quality and comprehensive reference map is provided as an open community resource for understanding human health and disease.",
  "groupsEdit": ["Admin"],
  "groupsView": ["Admin"],
  "info": {
    "Donors": 16,
    "Estimated_Cells": "528.1k",
    "Library_Construction_Method": "10X v2 sequencing",
    "Organ": "blood, immune system",
    "Organ_Part": "bone marrow, umbilical cord blood",
    "Paired_End": false,
    "Species": "Homo Sapian",
    "Specimens": 127
  },
  "name": "Census of Immune Cells",
  "s3path": "S3://dvc-portal-data/hca/8ca1f73e-b7e5-4b66-84b4-ca11ad14c528",
  "status": "success"
}
```

#### Configuration (Seed)

```JSON
{
  "id": "fa16a5e9-9c96-4ffb-b5ac-26eecc4b6eb7",
  "analysis": "d8c36499-4987-4ea4-9dc0-23c21ffde209",
  "desc": "Human Cell Atlas Inport",
  "groupsEdit": ["Admin"],
  "groupsView": ["Admin"],
  "info": {
    "date": 20200303,
    "source": "data.humancellatlas.org"
  },
  "name": "Human Cell Atlas Inport",
  "type": "imported"
}
```

#### Analysis (Seed)

```JSON
{
  "id": "d8c36499-4987-4ea4-9dc0-23c21ffde209",
  "analysis": ["9ccc42ac-e566-4f07-9ca2-135ac4ec2bb0"],
  "name": "HCA Loom File",
  "desc": "Analysis was preformed by the CZI HCA Data Portal To Create A Loom File",
  "files": ["data.loom"],
  "groupsView": ["Admin"],
  "groupsEdit": ["Admin"]
}
```

## Nextflow Script Execute + Modules

In production Nextflow pipelines will be triggered as a result of insert in the 'Configuration' dynamo db table when the type attribute equals 'one-time'.

By leveraging Nextflow DSL 2.0 we will create a generic method to read dataset, config, analysis objects as parameters in addition to numerically indexed input "parent" datasets
