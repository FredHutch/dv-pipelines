# Workflows


## Nextflow Script Execute + Modules

In production Nextflow pipelines will be triggered as a result of insert in the 'Configuration' dynamo db table when the type attribute equals 'one-time'.

By leveraging Nextflow DSL 2.0 we will create a generic method to read dataset, config, analysis objects as parameters in addition to numerically indexed input "parent" datasets
