## Run STARsolo workflow for processing 10X data
ml purge
ml STAR/2.7.3a-foss-2016b

## Base directory for analysis
BASE="/home/ramezqui/analysis/_TCRT"
cd $BASE

## Reference
REF="/home/ramezqui/reference/GENCODE/mouse"
REF_FASTA="${REF}/GRCm38.p6.genome.fa.gz"  # note it is genomic fasta not cDNA
REF_GTF="${REF}/gencode.vM24.chr_patch_hapl_scaff.annotation.gtf.gz"
REF_T2G="${REF}/gencode.vM24.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF

## Addon sequence
TCRT_FASTA="${BASE}/data-meta/tcrt-sequences_v3/addon.fa.gz"
TCRT_GTF="${BASE}/data-meta/tcrt-sequences_v3/addon.gtf.gz"

## Whitelist for technology
TENX_WHITELIST="${REF}/../../10x_whitelists/10xv2_whitelist.txt"

## Output dir for custom reference
REF_OUT="${BASE}/data-raw/STAR-ref_custom"
mkdir -p ${REF_OUT}

## Output dir for STAR results
STAR_OUT="${BASE}/data-raw/STAR-count"
mkdir -p ${STAR_OUT}

## Input dir of fastq files
FASTQ_GEX_DIR="${BASE}/data-raw/fastq_v2"


## Build custom reference index ------------------------------------------------
## Add together genomic fastas of chromosome sequences
zcat $TCRT_FASTA $REF_FASTA > ${REF_OUT}/custom.fa

## Add together GTF files of transcripts/genes/exons
zcat $REF_GTF | grep "\#" > ${REF_OUT}/custom.gtf
echo "##addendum: addition of custom TCRT gene 2020-04-29" >> ${REF_OUT}/custom.gtf
zcat $TCRT_GTF $REF_GTF | grep -v "\#" >> ${REF_OUT}/custom.gtf

## Standard transcriptome index
cd ${REF_OUT}

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
for DIR in $(ls -d ${FASTQ_GEX_DIR}/TCRT.41BB_spleen_*); do
    SAMPLE_ID=$(basename $DIR)

    echo "Aligning $SAMPLE_ID .."

    ## Writing an sbatch script
    echo '#!/bin/bash' > ${SAMPLE_ID}_star.sbatch
    echo "#SBATCH -p largenode" >> ${SAMPLE_ID}_star.sbatch
    echo "#SBATCH -n 1" >> ${SAMPLE_ID}_star.sbatch
    echo "#SBATCH -c 8" >> ${SAMPLE_ID}_star.sbatch
    echo "#SBATCH --mem=64000" >> ${SAMPLE_ID}_star.sbatch
    #SBATCH --mail-user=robert.amezquita@fredhutch.org
    #SBATCH --mail-type=ALL
    #SBATCH -D /fh/fast/gottardo_r/ramezqui_working/analysis/gottardo_hvtn097/data-raw
    #SBATCH -e /fh/fast/gottardo_r/ramezqui_working/slurm/master-%j.err
    #SBATCH -o /fh/fast/gottardo_r/ramezqui_working/slurm/master-%j.out
    echo "#SBATCH -J $SAMPLE_ID" >> ${SAMPLE_ID}_star.sbatch

    ## Switch directories
    echo "" >> ${SAMPLE_ID}_star.sbatch
    echo "mkdir -p ${STAR_OUT}/${SAMPLE_ID}" >> ${SAMPLE_ID}_star.sbatch
    echo "cd ${STAR_OUT}/${SAMPLE_ID}" >> ${SAMPLE_ID}_star.sbatch

    ## Standard + pre-mRNA + novel-splice junctions + velocity (see soloFeatures arg)
    ## - note that cDNA read must be first fastq input, 2nd is cell+UMI read
    ## - for filtering, CellRanger2.2 has 3 numbers:
    ##   - <expected cells> <max percentile for UMI count> <max to min ratio for UMI count>
    ## - the zcat could be replaced with "gunzip -c"
    ## - soloUMIfiltering and soloCBmatchWLtype are added to match more closely with CellRanger 3.x.x
    ##   UMI collapsing algorithm
    ## - bam file output has full list of standard tags
    echo "" >> ${SAMPLE_ID}_star.sbatch
    echo "STAR \
--runThreadN 8 \
--genomeDir ${REF_OUT} \
--readFilesCommand zcat \
--readFilesIn ${DIR}/*R2*.fastq.gz ${DIR}/*R1*.fastq.gz \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${TENX_WHITELIST} \
--soloCellFilter CellRanger2.2 10000 0.99 10 \
--soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloFeatures Gene GeneFull SJ Velocyto \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM" >> ${SAMPLE_ID}_star.sbatch

    chmod 755 ${SAMPLE_ID}_star.sbatch
    ./${SAMPLE_ID}_star.sbatch

    ## Alternately, if job queue isnt stacked..
    ## sbatch ${SAMPLE_ID}_star.sbatch
done