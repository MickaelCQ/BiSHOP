/*
    Reference: Picard MarkDuplicates via GATK4 (Broad Institute)
*/
process MARKDUPLICATES {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/aligned_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bai"), emit: bam_bai
    tuple val(meta), path("*.metrics.txt"),                         emit: metrics

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 4}g" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.markdup.bam \\
        -M ${meta.id}.markdup.metrics.txt \\
        --CREATE_INDEX true \\
        --REMOVE_DUPLICATES true
    """
}
