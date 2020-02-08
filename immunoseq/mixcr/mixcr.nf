fastq_files = Channel.fromPath(params.input.fastq_path + '/*.fq.gz')

process mixcr_analyze {
  echo false 

  publishDir "$params.output.folder/${fastq.getSimpleName()}"

  label 'gizmo_largenode'
  
  module 'MiXCR'
  
  scratch "/fh/scratch/delete30/warren_h/sravisha/MiXCR"

  input:
    path fastq from fastq_files
  
  output:
    path "${fastq.getSimpleName()}_mixcr.clonotypes.TRB.txt" into mixcr_clonotypes
    path "${fastq.getSimpleName()}.report" into mixcr_reports
    val task.exitStaus into count_status

  """
  mixcr analyze amplicon -s $params.analyze.species --starting-material $params.analyze.starting_material --5-end $params.analyze.five_end --3-end $params.analyze.three_end --adapters $params.analyze.adapters --receptor-type $params.analyze.receptor --report ${fastq.getSimpleName()}.report ${fastq} ${fastq.getSimpleName()}_mixcr
  """

}