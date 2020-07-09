# Aggregate Cellranger workflow

This workflow processes 2+ samples through cellranger count, and then runs cellranger aggr to aggregate the results



```
cd /Users/dnambi/Documents/GitHub/dv-pipelines/single-cell/cellranger-tds-irc/nf-simple 

nextflow run main.nf --source "s3://test-nextflow-data/sourceguid" --target "s3://test-nextflow-data/targetguid" --wfconfig config.json

```

Dataset 2 repro:

```

nextflow run main.nf --source "s3://test-nextflow-data/ade9b62f-ae3f-4817-81b4-ff53e1d73de1" --target "s3://test-nextflow-data/targetguid" --wfconfig dataset2.json

```


Local file setup:
```

nextflow run filelocal.nf --source "s3://test-nextflow-data/demo" --target "s3://test-nextflow-data/demo2" --wfconfig localfileconfig.json

```

S3 file setup:
```

nextflow run cloudfiles.nf --source "s3://test-nextflow-data/demo" --target "s3://test-nextflow-data/demo2" --wfconfig fileconfig.json

```

Cellranger with files staging

```

nextflow run cellrangercloud.nf --source "s3://test-nextflow-data/demo" --target "s3://test-nextflow-data/demo2" --wfconfig cellrangerfiles.json

```

































