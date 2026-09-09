/*
    Reference: Picard MarkDuplicates via GATK4 (Broad Institute)
    Rôle : Identification et élimination physique des doublons de PCR et des doublons optiques.
*/
process MARKDUPLICATES {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/aligned_bams" }, mode: 'link'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bai"), emit: bam_bai
    tuple val(meta), path("*.metrics.txt"),                         emit: metrics

    script:
    """
    # Options GATK MarkDuplicates :
    # -Xmx : Alloue la RAM Java en laissant 4 Go pour le système
    # --REMOVE_DUPLICATES true : Supprime physiquement les doublons PCR pour alléger le stockage
    # --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 : Distance en pixels adaptée aux patterned flowcells du NextSeq 2000
    # --CREATE_INDEX true : Crée automatiquement l'index .bai
    # --VALIDATION_STRINGENCY SILENT : Ignore les avertissements bénins
    gatk --java-options "-Xmx${task.memory.toGiga() - 4}g -XX:+UseParallelGC" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.markdup.bam \\
        -M ${meta.id}.markdup.metrics.txt \\
        --REMOVE_DUPLICATES true \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT
    """
}
