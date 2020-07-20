#!/usr/bin/env nextflow
x_val = "sample2"
input_file_loc = "/root/doubletest/sourcefile.txt"

values = Channel.of( [x_val, input_file_loc] )

process TEST_MAP {
  echo true
  publishDir '/root/doubletest/outputtest', mode: 'copy'

  input:
    tuple val(sample), path('source') from values

  output:
    tuple val(sample), path("${sample}/*") into pub_ch

  script:
    """
    mkdir -p $sample

    touch $sample/testfile.txt
    echo "line1" >> $sample/testfile.txt
    echo 'ignore this' >> $sample/testfile.txt
    cat $source >> $sample/testfile.txt
    echo "Input file contents:"
    cat $source
    ls
    """
}
