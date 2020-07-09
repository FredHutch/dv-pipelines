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

## Input dir of fastq files
#FASTQ_GEX_DIR="s3://test-nextflow-data/40d55f14-5ae8-4e7a-af9b-85d1122a0aaa"
FASTQ_GEX_DIR="${BASE}/fastq"

## See https://github.com/nf-core/scrnaseq/blob/master/docs/usage.md

REF="${BASE}/reference/GENCODE/mouse"
#REF="s3://hutch.cloud.dev/reference-data/GENCODE/genomes/mouse"
REF_FASTA="${REF}/GRCm38.p6.genome.fa.gz"  # note it is genomic fasta not cDNA
REF_GTF="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz"
gunzip $REF_FASTA
gunzip $REF_GTF


REF="${BASE}/reference/GENCODE/mouse"
REF_FASTA="${REF}/GRCm38.p6.genome.fa"
REF_GTF="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf"
REF_T2G="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF

nextflow run nf-core/scrnaseq -r 1.0.0 --reads "${FASTQ_GEX_DIR}/*_R{1,2}.fastq.gz" \
 -profile docker --aligner 'star' --fasta $REF_FASTA --gtf $REF_GTF
