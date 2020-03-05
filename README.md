# Nextflow pipelines for scRNA-Seq and ImmunoSeq analysis

The repository contains nextflow scripts to run different workflows for scRNA-Seq analysis as well as MiXCR workflow for AIRR-Seq data analysis. The scRNA-Seq workflows are divide into five logical steps

1. Count matrix generation
2. Preprocessing
3. Clustering
4. Cell-type classification
5. Differential expression analysis
6. Trajectory analysis


## scRNA-Seq and AIRR-Seq tools incorporated

1. MiXCR
2. Monocle3
3. Seurat
4. Scanpy
5. Scater
6. Scran
7. DropletUtils

# Contributing Workflows

Contributing workflows requires the following files.

## EXAMPLE.FORM.*

Specifics on all files that start with TEST.FORM can be be found here:

Documentation Site: <https://react-jsonschema-form.readthedocs.io/en/latest/advanced-customization/>

Online Test Harness: <https://cybertec-postgresql.github.io/rjsf-material-ui/>

## EXAMPLE.FORM.JSON

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

## EXAMPLE.FORM.UI.JSON

This file can used to further customize the form.  It can include help text, enumeration values, choice of UI Components
If you are comfortable with a basic form then this file can simply include an empty JSON object (EG {})

```JSON
{
  "n_comps": {
    "ui:widget": "updown"
  }
}
```

## EXAMPLE.FORM.VALUES.JSON

This file contains a json object with the default values to render in the form.
When the user specifies parameters to run they will be serialized in the same format + placed in the nextflow working directory as
"params.json"

```JSON
{
  "zero_center": false,
  "n_comps": 50
}
```



## Dynamo Records To Run Analysis

Analysis (Test Dataset)

```JSON
{
  "id": "9ccc42ac-e566-4f07-9ca2-135ac4ec2bb0",
  "analysis": [],
  "name": "Test Analysis",
  "desc": "Basic Loom Test Analysis",
  "info": {},
  "runtime": "NEXTFLOW",
  "executable": "https://github.com/FredHutch/dv-pipelines/test.nf",
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

Analysis (Source Dataset)

```JSON
{
  "id": "d8c36499-4987-4ea4-9dc0-23c21ffde209",
  "analysis": ["9ccc42ac-e566-4f07-9ca2-135ac4ec2bb0"],
  "name": "HCA Loom File",
  "desc": "Analysis was preformed by the CZI HCA Data Portal To Create A Loom File",
  "info": "{}",
  "runtime": "NA",
  "executable": "NA",
  "schema": "{}",
  "files": ["data.loom"],
  "groupsView": ["Admin"],
  "groupsEdit": ["Admin"]

}
```

Configuration  (Source Dataset)

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
  "parameters": "{}",
  "type": "IMPORTED"
}
```

Dataset (Source Dataset)

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

