texts = Channel.from("AWS", "Nextflow")

process hello {
    // directives
    // a container images is required
    container "ubuntu:latest"

    // compute resources for the Batch Job
    cpus 1
    memory '512 MB'

    input:
    val text from texts

    output:
    file 'hello2.txt'

    """
    echo "Hello $text" > hello.txt
    """
}