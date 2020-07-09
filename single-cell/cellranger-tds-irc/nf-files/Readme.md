# NF File-handling experiments

Create files in a folder, get them output to a different process

```
cd /Users/dnambi/Documents/GitHub/dv-pipelines/single-cell/cellranger-tds-irc/nf-files/1-copyoutputs/ 

nextflow run copyoutputs.nf --source "s3://test-nextflow-data/sourceguid" --target "s3://test-nextflow-data/targetguid" --wfconfig localfileconfig.json

```

Do a publish test

```
cd /Users/dnambi/Documents/GitHub/dv-pipelines/single-cell/cellranger-tds-irc/nf-files/2-publish/

nextflow run publish.nf --source "s3://test-nextflow-data/sourceguid" --target "s3://test-nextflow-data/scratch/dnambitests/2-publish" --wfconfig config.json

```


```
ssh -i "dnambi-vm-key-sttr.pem" ubuntu@ec2-54-212-60-197.us-west-2.compute.amazonaws.com


```



















