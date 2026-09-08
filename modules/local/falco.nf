/*
    Reference: Falco - C++ port of FastQC (https://github.com/gt1/falco)
*/
process FALCO {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/falco/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt") , emit: txt

    script:
    """
    falco --threads ${task.cpus} ${reads[0]} ${reads[1]}
    """
}
