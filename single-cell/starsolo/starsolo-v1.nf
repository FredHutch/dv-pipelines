#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

// "${BASE}/data-raw/STAR-count"
target_path = "$params.s3bucket" + "/"
target_path_val = Channel.value("$params.s3bucket" + "/")


// "/home/ramezqui/analysis/_BCMA"
source_path = wfi.parents[0].s3path + "/"
source_path_val = source_path

//Input parameters
/// Reference data

// "/home/ramezqui/reference/GENCODE/human"
reference_path = wfi.parameters.input.reference_path 


// "${REF}/GRCh38.p13.genome.fa.gz"  # note it is genomic fasta not cDNA
ref_fasta = wfi.parameters.input.reference_path + wfi.parameters.input.reference_fasta


// "${REF}/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
ref_gtf = wfi.parameters.input.reference_path + wfi.parameters.input.reference_gtf


// "${REF}/gencode.v33.chr_patch_hapl_scaff.annotation.t2g" # parsed from GTF
ref_t2g = wfi.parameters.input.reference_path + wfi.parameters.input.reference_t2g

// Addon sequence

//"${BASE}/data-meta/cart-sequence/addon.fa.gz" 
cart_fasta = wfi.parameters.input.cart_fasta

// "${BASE}/data-meta/cart-sequence/addon.gtf.gz"
cart_gtf = wfi.parameters.input.cart_gtf


// "${BASE}/data-meta/cart-sequence/tx2gene.txt"
cart_t2g = wfi.parameters.input.cart_t2g


// Whitelist for technology
// "${REF}/../../10x_whitelists/10xv2_whitelist.txt"
tenx_whitelist = wfi.parameters.input.tenx_whitelist

// Input dir of fastq files
// "${BASE}/data-raw/fastq/samples_GEX-only"
gex_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)




// *maybe* necessary
// "${BASE}/data-raw/STAR-ref_custom"
output_ref = target_path + "STAR-ref-custom/"
sjdbOverhang = wfi.parameters.input.sjdbOverhang


// Flags
run_create_ref_index = Channel.from(wfi.parameters.input.create-ref-index)
run_genome_create = Channel.from(wfi.parameters.input.create-genome)



process CREATE_CUSTOM_REF_INDEX {
  echo true
  publishDir ?

  input:
    cart_fasta
    ref_fasta
    output_ref
    ref_gtf
    cart_gtf

  output:
    custom.fa
    custom.gtf

  when:
    run_create_ref_index == 1

  script:
    """
    // Build custom reference index ------------------------------------------------
    // Add together genomic fastas of chromosome sequences
    zcat $CART_FASTA $REF_FASTA > custom.fa
    // gzip custom.fa

    // Add together GTF files of transcripts/genes/exons
    zcat $REF_GTF | grep "\#" > custom.gtf
    echo "//addendum: addition of custom CART gene targeting BCMA protein 2020-04-29" >> custom.gtf
    zcat $CART_GTF $REF_GTF | grep -v "\#" >> custom.gtf
    """
}



process CREATE_GENOME {
  echo true
  publishDir ?

  input:
    custom.fa
    custom.gtf

  output:
    work/.* from work into sample_genomes // split by folder

  when:
    run_genome_create == 1

  script:
    """
    mkdir -p work
    // Standard transcriptome index
    STAR --runMode genomeGenerate \
         --runThreadN 8 \
         --genomeDir work \
         --genomeFastaFiles custom.fa \
         --sjdbGTFfile /custom.gtf \
         --sjdbOverhang {sjdbOverhang}
    """
}






// Standard + pre-mRNA + novel-splice junctions + velocity (see soloFeatures arg)
// - note that cDNA read must be first fastq input, 2nd is cell+UMI read
// - for filtering, CellRanger2.2 has 3 numbers:
//   - <expected cells> <max percentile for UMI count> <max to min ratio for UMI count>
// - the zcat could be replaced with "gunzip -c"
// - soloUMIfiltering and soloCBmatchWLtype are added to match more closely with CellRanger 3.x.x
//   UMI collapsing algorithm
// - bam file output has full list of standard tags
// needs 8 cores and 64GB RAM
process QUANTIFY_SAMPLES {
  echo true
  publishDir ?

  input:
    val sample_dir from gex_fastq_files
    val sample_id from sample_genomes
    val tenx_whitelist
    val expected_cells

  output:
    work/.*

  script:
    """
    mkdir -p work
    ## Perform quantification per sample -------------------------------------------
    SAMPLE_ID=\$(basename ${sample_dir})
    echo "Aligning \$SAMPLE_ID .."
    STAR \
    --runThreadN {CPU_THREADS} \
    --genomeDir work \
    --readFilesCommand zcat \
    --readFilesIn ${DIR}/*R2*.fastq.gz ${DIR}/*R1*.fastq.gz \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ${tenx_whitelist} \
    --soloCellFilter CellRanger2.2 {expected_cells} 0.99 10 \
    --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
    """
}

