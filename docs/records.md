# Records


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
