## Run STARsolo workflow for processing 10X data

# PART 1: SET UP THE DATA AND BASE VARIABLES
## Base directory for analysis
BASE="/scratch/workflow"
cd $BASE


nextflow run nf-core/scrnaseq -r 1.0.0 -c custom.config -profile docker
