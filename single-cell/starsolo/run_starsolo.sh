## Run STARsolo workflow for processing 10X data

ml STAR/2.7.3a-foss-2016b

## Base directory for analysis
BASE="/home/ramezqui/analysis/_BCMA"
cd $BASE

## Reference 
REF="/home/ramezqui/reference/GENCODE/human"
REF_FASTA="${REF}/GRCh38.p13.genome.fa.gz"  # note it is genomic fasta not cDNA
  ## DEVN: what's the difference between genomic fasta and cDNA?
REF_GTF="${REF}/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
REF_T2G="${REF}/gencode.v33.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF
## DEVN: chr_patch means what? Hapl scaff? 

## Addon sequence
CART_FASTA="${BASE}/data-meta/cart-sequence/addon.fa.gz"
CART_GTF="${BASE}/data-meta/cart-sequence/addon.gtf.gz"
CART_T2G="${BASE}/data-meta/cart-sequence/tx2gene.txt"

## Whitelist for technology
## DEVN: what is the whitelist for??
TENX_WHITELIST="${REF}/../../10x_whitelists/10xv2_whitelist.txt"

## Output dir for custom reference
REF_OUT="${BASE}/data-raw/STAR-ref_custom"
mkdir -p ${REF_OUT}
cd ${REF_OUT}

## Output dir for STAR results
STAR_OUT="${BASE}/data-raw/STAR-count"
mkdir -p ${STAR_OUT}

## Input dir of fastq files
FASTQ_GEX_DIR="${BASE}/data-raw/fastq/samples_GEX-only"


## Build custom reference index ------------------------------------------------
## Add together genomic fastas of chromosome sequences
zcat $CART_FASTA $REF_FASTA > ${REF_OUT}/custom.fa
## gzip ${REF_OUT}/custom.fa

## Add together GTF files of transcripts/genes/exons
zcat $REF_GTF | grep "\#" > ${REF_OUT}/custom.gtf
echo "##addendum: addition of custom CART gene targeting BCMA protein 2020-04-29" >> ${REF_OUT}/custom.gtf
zcat $CART_GTF $REF_GTF | grep -v "\#" >> ${REF_OUT}/custom.gtf

## Standard transcriptome index
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${REF_OUT} \
     --genomeFastaFiles ${REF_OUT}/custom.fa \
     --sjdbGTFfile ${REF_OUT}/custom.gtf \
     --sjdbOverhang 120
     ## DEVN: why the sjdbOverhang of 120? 

## Perform quantification per sample -------------------------------------------
## Run sbatch jobs in a runs directory to save slurm output
mkdir -p ${BASE}/runs
cd ${BASE}/runs

for DIR in $(ls -d ${FASTQ_GEX_DIR}/*); do
    SAMPLE_ID=$(basename $DIR)

    echo "Aligning $SAMPLE_ID .."

    ## Writing an sbatch script
    echo '#!/bin/bash' > ${SAMPLE_ID}.sbatch
    echo "#SBATCH -p largenode" >> ${SAMPLE_ID}.sbatch
    echo "#SBATCH -n 1" >> ${SAMPLE_ID}.sbatch
    echo "#SBATCH -c 8" >> ${SAMPLE_ID}.sbatch
    echo "#SBATCH --mem=64000" >> ${SAMPLE_ID}.sbatch
    #SBATCH --mail-user=robert.amezquita@fredhutch.org
    #SBATCH --mail-type=ALL
    #SBATCH -D /fh/fast/gottardo_r/ramezqui_working/analysis/gottardo_hvtn097/data-raw
    #SBATCH -e /fh/fast/gottardo_r/ramezqui_working/slurm/master-%j.err
    #SBATCH -o /fh/fast/gottardo_r/ramezqui_working/slurm/master-%j.out
    echo "#SBATCH -J $SAMPLE_ID" >> ${SAMPLE_ID}.sbatch

    ## Switch directories
    echo "" >> ${SAMPLE_ID}.sbatch
    echo "mkdir -p ${STAR_OUT}/${SAMPLE_ID}" >> ${SAMPLE_ID}.sbatch
    echo "cd ${STAR_OUT}/${SAMPLE_ID}" >> ${SAMPLE_ID}.sbatch

    ## Standard + pre-mRNA + novel-splice junctions + velocity (see soloFeatures arg)
    ## - note that cDNA read must be first fastq input, 2nd is cell+UMI read
    ## - for filtering, CellRanger2.2 has 3 numbers:
    ##   - <expected cells> <max percentile for UMI count> <max to min ratio for UMI count>
    ## - the zcat could be replaced with "gunzip -c"
    ## - soloUMIfiltering and soloCBmatchWLtype are added to match more closely with CellRanger 3.x.x
    ##   UMI collapsing algorithm
    ## - bam file output has full list of standard tags
    echo "" >> ${SAMPLE_ID}.sbatch
    echo "STAR \
--runThreadN 8 \
--genomeDir ${REF_OUT} \
--readFilesCommand zcat \
--readFilesIn ${DIR}/*R2*.fastq.gz ${DIR}/*R1*.fastq.gz \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${TENX_WHITELIST} \
--soloCellFilter CellRanger2.2 10000 0.99 10 \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--soloFeatures Gene GeneFull SJ Velocyto \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM" >> ${SAMPLE_ID}.sbatch
# DEVN: So this does count quantification? 
# DEVN: soloType CB_UMI_Simple means what? 
# --soloCellFilter - 
# one of the key problem with scRNAseq is you don't know how many
# cells end up getting sequenced. Typically they expect 10K cells max
# Need to figure out how many cells where actually sequenced
# STAR always shows the raw matrix, cell barcodes + UMIs + read itself
# it's split by cell barcode, generally maps to a single droplet
# helps to figure out what the valid cells is
# number to get from the experimentalist is est. # of cells from that sample
# the default seems to be 3K, 0.99, 10. 
# INPUT: just ask for the number of expected cells from the experimentalist
# 

    chmod 755 ${SAMPLE_ID}.sbatch
##    ./${SAMPLE_ID}.sbatch

    ## Alternately, if job queue isnt stacked..
    ## sbatch ${SAMPLE_ID}.sbatch
done

