#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)
source_path = "$params.s3bucket" + "/"

//Input parameters
/// Reference data
reference_path = wfi.parameters.reference_path
gex_reference = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)

/// Workflow data
metadata = Channel.fromPath(source_path + wfi.parameters.input.metadata)
gex_fastq_paths = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
vdj_fastq_paths = Channel.fromPath(source_path + wfi.parameters.input.vdj_fastq_path)
study_id = Channel.from(wfi.parameters.input.study_id)
modes = Channel.from(wfi.parameters.aggr.modes)
run_gex = Channel.from(wfi.parameters.input.gex)
run_vdj = Channel.from(wfi.parameters.input.vdj)
fastq_type = Channel.from(wfi.parameters.count.fastq_type)

//Split input channels as needed
run_gex.into{map_gex ; count_gex ; aggr_gex ; mat_gex}
metadata.into{gex_metadata ; vdj_metadata}
run_vdj.into{map_vdj ; ana_vdj}
study_id.into{gex_study_id; vdj_study_id}




process PRINT_ALL {
  echo true
  input:
    file gex_meta_file from gex_metadata
    each mode from modes
    val gex_fastq_path from gex_fastq_paths
    val study from gex_study_id
    val map from map_gex
    val vdj_fastq_path from vdj_fastq_paths
    val analyze from ana_vdj
    val gex_reference from gex_reference
    val vdj_reference from vdj_reference
    val fastq_type from fastq_type

  script:
    """
    echo "Mode $mode"
    echo "gex_fastq_path $gex_fastq_path"
    echo "study $study"
    echo "map $map"
    echo "vdj_fastq_path $vdj_fastq_path"
    echo "analyze $analyze"
    echo "gex_reference $gex_reference"
    echo "vdj_reference $vdj_reference"
    echo "fastq_type $fastq_type"
    cat $gex_meta_file
    """

}
