## Run STARsolo workflow for processing 10X data

# PART 1: SET UP THE DATA AND BASE VARIABLES
## Base directory for analysis
BASE="/scratch/workflow"
cd $BASE

## Input dir of fastq files
#FASTQ_GEX_DIR="s3://test-nextflow-data/40d55f14-5ae8-4e7a-af9b-85d1122a0aaa"
FASTQ_GEX_DIR="${BASE}/fastq"


REF="${BASE}/reference/GENCODE/mouse"
REF_FASTA="${REF}/GRCm38.p6.genome.fa"
REF_GTF="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf"
REF_T2G="${REF}/gencode.vM25.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF

nextflow run nf-core/scrnaseq -r 1.0.0 --reads "${FASTQ_GEX_DIR}/*_R{1,2}.fastq.gz" \
 -profile docker --aligner 'star' --fasta $REF_FASTA --gtf $REF_GTF
