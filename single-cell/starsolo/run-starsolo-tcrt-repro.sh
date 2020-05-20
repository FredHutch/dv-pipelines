## Run STARsolo workflow for processing 10X data

# PART 1: SET UP THE DATA AND BASE VARIABLES
mkdir -p /scratch/workflow

## Base directory for analysis
BASE="/scratch/workflow"
cd $BASE

## Get the reference data
aws s3 cp s3://hutch.cloud.dev/reference-data/GENCODE/genomes/mouse/ reference/GENCODE/mouse/ --recursive
aws s3 cp s3://hutch.cloud.dev/reference-data/10x_whitelists/ reference/10x_whitelists/ --recursive

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

## Whitelist for technology
TENX_WHITELIST="${BASE}/reference/10x_whitelists/10xv2_whitelist.txt"

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



STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${REF_OUT} \
     --genomeFastaFiles ${REF_OUT}/custom.fa \
     --sjdbGTFfile ${REF_OUT}/custom.gtf \
     --sjdbOverhang 120



## Perform quantification per sample -------------------------------------------
## Run sbatch jobs in a runs directory to save slurm output
mkdir -p ${BASE}/runs
cd ${BASE}/runs

##for DIR in $(ls -d ${FASTQ_GEX_DIR}/*); do
for DIR in $(ls ${FASTQ_GEX_DIR}/*.fastq.gz | grep R1 | sed 's/_R1.fastq.gz//g'); do
    SAMPLE_ID=$(basename $DIR)

    echo "Aligning $SAMPLE_ID .."

    ## Switch directories
    mkdir -p ${STAR_OUT}/${SAMPLE_ID}
    cd ${STAR_OUT}/${SAMPLE_ID}

    ## Standard + pre-mRNA + novel-splice junctions + velocity (see soloFeatures arg)
    ## - note that cDNA read must be first fastq input, 2nd is cell+UMI read
    ## - for filtering, CellRanger2.2 has 3 numbers:
    ##   - <expected cells> <max percentile for UMI count> <max to min ratio for UMI count>
    ## - the zcat could be replaced with "gunzip -c"
    ## - soloUMIfiltering and soloCBmatchWLtype are added to match more closely with CellRanger 3.x.x
    ##   UMI collapsing algorithm
    ## - bam file output has full list of standard tags
    STAR \
    --runThreadN 8 \
    --genomeDir ${REF_OUT} \
    --readFilesCommand zcat \
    --readFilesIn ${DIR}*R2*.fastq.gz ${DIR}*R1*.fastq.gz \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ${TENX_WHITELIST} \
    --soloCellFilter CellRanger2.2 10000 0.99 10 \
    --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
done