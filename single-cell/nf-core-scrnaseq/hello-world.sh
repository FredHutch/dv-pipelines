## Run STARsolo workflow for processing 10X data

# PART 1: SET UP THE DATA AND BASE VARIABLES
mkdir -p /scratch/workflow

## Base directory for analysis
BASE="/scratch/workflow"
cd $BASE

## Get the reference data
aws s3 cp s3://hutch.cloud.dev/reference-data/GENCODE/genomes/mouse/ reference/GENCODE/mouse/ --recursive

## Get the rest of the data
aws s3 cp s3://test-nextflow-data/40d55f14-5ae8-4e7a-af9b-85d1122a0aaa/ /scratch/workflow/ --recursive


## Reference
REF="${BASE}/reference/GENCODE/mouse"
REF_FASTA="${REF}/GRCm38.p6.genome.fa.gz"  # note it is genomic fasta not cDNA
REF_GTF="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz"
REF_T2G="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF

## Addon sequence
TCRT_FASTA="${BASE}/addon/addon.fa"
TCRT_GTF="${BASE}/addon/addon.gtf"
gzip -f $TCRT_FASTA
gzip -f $TCRT_GTF
TCRT_FASTA="${BASE}/addon/addon.fa.gz"
TCRT_GTF="${BASE}/addon/addon.gtf.gz"

## Output dir for custom reference
REF_OUT="${BASE}/data-raw/STAR-ref_custom"
mkdir -p ${REF_OUT}

## Output dir for STAR results
STAR_OUT="${BASE}/data-raw/STAR-count"
mkdir -p ${STAR_OUT}

## Input dir of fastq files
FASTQ_GEX_DIR="${BASE}/fastq"



## Build custom reference index ------------------------------------------------
## Add together genomic fastas of chromosome sequences
zcat $TCRT_FASTA $REF_FASTA > ${REF_OUT}/custom.fa

## Add together GTF files of transcripts/genes/exons
zcat $REF_GTF | grep "\#" > ${REF_OUT}/custom.gtf
echo "##addendum: addition of custom TCRT gene 2020-04-29" >> ${REF_OUT}/custom.gtf
zcat $TCRT_GTF $REF_GTF | grep -v "\#" >> ${REF_OUT}/custom.gtf


## See https://github.com/nf-core/scrnaseq/blob/master/docs/usage.md


nextflow run nf-core/scrnaseq -r 1.0.0 --reads "${FASTQ_GEX_DIR}/*_R{1,2}.fastq.gz" \
   --fasta ${REF_OUT}/custom.fa --gtf ${REF_OUT}/custom.gtf -profile docker --aligner STAR
