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
reference_path = wfi.parameters.input.reference_path
gex_reference = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)

gex_reference2 = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference2 = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)


/// Workflow data
metadata = Channel.fromPath(target_path + wfi.parameters.input.metadata)
gex_fastq_paths = Channel.of(source_path + wfi.parameters.input.gex_fastq_path)
gex_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
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



process DEBUG {
  echo true
  scratch = '/opt/work'

  input: 
    val gex_reference2
    file gex_ref from gex_reference
    val vdj_reference2

  script:
    """
    echo "GEX Reference is $gex_reference2"
    
    mkdir ref
    tar -zxvf $gex_ref -C ref
    du -sh ref
    REFPATH="\$(ls ref/)"
    ls -lh ref/\$REFPATH
    """
}