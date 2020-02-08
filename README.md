## Nextflow: Getting started 

1. Load nextflow on gizmo

```
ml nextflow
```

2. Hello world

Lets now try the hello world script from the next getting started documentation. The workflow defines two processes, splitLetters and convertToUpper. The parameters supplied to the workflow are params.str

```
#!/usr/bin/env nextflow

params.str = "Hello world!"

process splitLetters {
	
	output:
	file 'chunk_*' into letters

	"""
	printf '${params.str}' | split -b 6 - chunk_
	"""

}

process convertToUpper {
	
	input:
	file x from letters.flatten()

	output:
	stdout result

	"""
	cat $x | tr '[a-z]' '[A-Z]'
	"""
}

results.view { it.trim() }
```


Run this script using the command

```
nextflow run tutorial.nf
```
