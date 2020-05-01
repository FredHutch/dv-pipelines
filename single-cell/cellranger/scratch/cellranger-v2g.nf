#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)
target_path = "$params.s3bucket" + "/"
target_path_val = Channel.value("$params.s3bucket" + "/")
source_path = wfi.parents[0].s3path + "/"
source_path_val = source_path


//Input parameters
/// Reference data
reference_path = wfi.parameters.reference_path
gex_reference = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)

/// Workflow data
metadata = Channel.fromPath(target_path + wfi.parameters.input.metadata)
gex_fastq_paths = Channel.of(source_path + wfi.parameters.input.gex_fastq_path)
gex_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
gex_fastq_files2 = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
vdj_fastq_paths = Channel.of(source_path + wfi.parameters.input.vdj_fastq_path)
vdj_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.vdj_fastq_path)
study_id = Channel.from(wfi.parameters.input.study_id)
modes = Channel.from(wfi.parameters.aggr.modes)
run_gex = Channel.from(wfi.parameters.input.gex)
run_vdj = Channel.from(wfi.parameters.input.vdj)
fastq_type = Channel.value(wfi.parameters.count.fastq_type)

//Split input channels as needed
run_gex.into{map_gex ; count_gex ; aggr_gex ; mat_gex}
metadata.into{gex_metadata ; vdj_metadata}
run_vdj.into{map_vdj ; ana_vdj}
study_id.into{gex_study_id; vdj_study_id}


process TENX_COUNT {
  echo true 
  scratch = '/opt/work'

  input:
    file sample_file from gex_fastq_files.collect()
    val filepaths from gex_fastq_files2.collect()
    val run_count from count_gex
    val target_path_val
    val source_path_val
  
  output:
    val task.exitStaus into count_status

  when:
    run_count == 1

  script:
    """
    echo "Source folder is is $source_path_val"
    echo "Sample files are $filepaths"

    mkdir fastq
    mv $sample_file fastq
    du -sh fastq
    ls -lh fastq | wc -l > filesize.log
    ls -lh fastq >> filesize.log
    cat filesize.log
    """
}



