#!/usr/bin/env nextflow
//Input parameters
target_path = "$params.s3target"
source_path = "$params.s3source"

input_files = Channel.fromPath( source_path )

process WORD_COUNT {
  echo true
  publishDir target_path, mode: 'copy'

  input:
    path x from input_files

  output:
    file "filelist.txt" into hdf5_ch

  script:
    """
    echo "Local files are:"
    ls -R * > filelist.txt
    cat filelist.txt
    """
}
