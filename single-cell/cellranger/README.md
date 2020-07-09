# Cellranger Nextflow workflow

This cellranger Nextflow workflow is part of the Single cell RNASeq workflow developed for the Warren Lab. The following document outlines the input and output specifications, data formatting requirements and runtime enviroment details for the workflow. 


## Execution

```
nextflow run cellranger.nf --wfconfig 'config.json'
```

## Inputs:

All the inputs for the analysis are read out of a config file, the structure of the config file is listed below.

```{nextflow}
params {
    input {
      metadata = 'ZhengSorted_metadata.csv'
      gex_reference = 'refdata-cellranger-hg19-3.0.0'
      vdj_reference = 'refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0'
      fastq_paths = 'ZhengSorted_10X/Fastq'
      study_id = 'ZhengSorted_10X'
      gex = true
      vdj = false
    }

    output {
      folder = "ZhengSorted_10X/Results"
    }

    count {
      fastq_type = 'demux' //['mkfastq', 'demux', 'bcl2fastq']
    }

    aggr {
      modes = "mapped" //['mapped', 'none']
    }
}
```

The parameters file is divided into four sections.

1. Input

   The input section contains the path to the metadata file, GEX and VDJ reference genome paths, path to the folders containing Fastq files, study name and flags to determine if VDJ and GEX analysis needs to be done.

   i. Metadata

      The metadata file outlines all the library information for the given study and must contain the columns shown in the example below. Any additional columns listed in the file will be added to the SingleCellExperiment object column data in later workflows. 

    | repoName | sampleName | patientID | expected_cells | chemistry | nucliecAcid | locus | platform | indices | sortingCT | tissueType | VDJType
    -|----------|------------|-----------|----------------|-----------|-------------|-------|----------|---------|-----------|------------|--------
    | S1_Monocytes | Monocytes | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | A1 | CD14pMonocytes | PBMC | NA
    S1_Bcells | Bcells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | A2 | CD19PBCells | PBMC | NA
    S1_Progenitor | Progenitor | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | E1 | CD34pCells | PBMC | NA
    S1_HelperTCells | HelperTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | E2 | CD4pHelperTCells | PBMC | NA
    S1_RegulatoryTCells | RegulatoryTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | D1 | CD4p_CD25pRegulatoryCells | PBMC | NA
    S1_NaiveTCells | NaiveTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | D2 | CD4p_CD45RAp_CD25nNaiveTCells | PBMC | NA
    S1_MemoryTCells | MemoryTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | B1 | CD4p_CD45ROpMemoryTCells | PBMC | NA
    S1_NKCells | NKCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | B2 | CD56pNKCells | PBMC | NA
    S1_CytotoxicTCells | CytotoxicTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | G1 | CD8pCytotoxicTCells | PBMC | NA
    S1_NaiveCytotoxicTCells | NaiveCytotoxicTCells | S1 | 3000 | threeprime | cDNA | 3primeGEX | 10XGenomics | G2 | CD8p_CD45RApNaiveCytotoxicTCells | PBMC | NA

    The VDJType column would either be `T cell` or `B cell` if VDJ sequencing was also done on these sample, if not the column value will be `NA`.

    ii. GEX and VDJ reference files

    The reference genome files for GEX and VDJ analysis can be downloaded from the 10X Genomics website using the links mentioned below:

    ```
    GEX reference: http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
    VDJ reference: http://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0.tar.gz
    ```

    iii. Path to folder containing FASTQ files

    The folder structure of the FASTQ inputs can vary depending on the method used to generate the FASTQ from the raw sequencing files. The current pipeline FASTQ files generated using two methods `demux` and `mkfastq`. The recommended folder structure is mentioned below can be found at [this link.](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input) The current version of the pipeline was developed considering cellranger `mkfastq` outputs, and partially tailored to accept cellranger `demux` output. Other options will be made avaiable in future versions.

2. Output
   
   All results from the cellranger Nextflow will be stored within subfolders under this output path.

3. Count

   Currently the only parameter specified here is the method used to generate the FASTQ files, the additional run details such as sequencing chemistry and expected number of cells are obtained from the metadata file. Runtime characteristics such as memory and number of CPUs are obtain from the process information specified in the nextflow.config file. 

4. Aggregate

   Specifies what type of normalization should be performed during the aggregation step.
    

## Outputs:

The pipeline generates three folders in the output folder:

1. Counts

   Contains the count matrices for each samples in the analysis. This serves as the input for most of the other nextflow scripts for downstream analysis

2. Metadata
   
   The metadata folder contain the samplesheets for the VDJ and GEX analysisas well as additional information such files paths for count matrices and VDJ clonotypes and contig files. These sample metadta files serve as the input to the preprocessing nextflow that generate the Single Cell Experiment objects for each sample. All information from the sample metadata files is added to the colData of the Single Cell Experimentß

3. VDJ
   
   The VDJ folder contains the results from the cellranger VDJ pipeline when available. The folder contains the clonotypes and contig. When avaiable the VDJ data is incorporated into the Single Cell Experiment object in the preprocessing nextflow

