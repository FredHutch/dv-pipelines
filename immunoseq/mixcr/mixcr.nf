fastq_files = Channel.fromPath(params.input.fastq_path + '/*.fq.gz')
minsumscores = Channel.from(params.analyze.minsumscore)
maxhits = Channel.from(params.analyze.maxhits)
minseqlength = Channel.from(params.analyze.minseqlength)

process mixcr_analyze {
  echo false 

  publishDir "$params.output.folder/${fastq.getSimpleName()}", mode: 'copy'

  label 'gizmo_largenode'
  
  module 'MiXCR'
  
  scratch "/fh/scratch/delete30/warren_h/sravisha/MiXCR"

  input:
    val fastq from fastq_files
    each minSumScore from minsumscores
    each maxHits from maxhits
    each minSeqLength from minseqlength
  
  output:
    path "${fastq.getSimpleName()}_mixcr.clonotypes.TRB.txt" into mixcr_clonotypes
    path "${fastq.getSimpleName()}.report" into mixcr_reports
    val task.exitStaus into count_status

  """
  mixcr analyze amplicon -s $params.analyze.species --starting-material $params.analyze.starting_material --5-end $params.analyze.five_end \
                        --3-end $params.analyze.three_end --adapters $params.analyze.adapters --receptor-type $params.analyze.receptor \
                        --align "-OminSumScore=${minSumScore} -OmaxHits=${maxHits}"  --assemble "-OminimalClonalSequenceLength=${minSeqLength}" \
                        --report ${fastq.getSimpleName()}.report ${fastq} ${fastq.getSimpleName()}_mixcr 
  """

}