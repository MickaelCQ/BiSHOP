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


/*
    Reference: Picard MarkDuplicates via GATK4 (Broad Institute)
    Rôle : Identification et élimination physique des doublons de PCR et des doublons optiques.
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
    gatk --java-options "-Xmx${task.memory.toGiga() - 4}g -XX:+UseParallelGC" MarkDuplicates \\
        # -Xmx : Réserve la RAM Java allouée par SLURM en laissant 4 Go pour le système hôte.
        -I ${bam} \\
        -O ${meta.id}.markdup.bam \\
        -M ${meta.id}.markdup.metrics.txt \\
        # Fichier texte de métriques contenant le taux de duplication exact analysé par MultiQC.
        --REMOVE_DUPLICATES true \\
        # Supprime physiquement les reads dupliqués du BAM pour alléger le stockage NFS et garantir qu'aucun
        # caller (notamment FreeBayes et les outils CNV) ne prenne un clone PCR pour un vrai allèle.
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        # Distance en pixels adaptée aux flowcells (Patterned Flowcells) du NextSeq 2000 et NovaSeq
        # pour différencier les doublons optiques de surface des vrais fragments indépendants.
        --CREATE_INDEX true \\
        # Génère automatiquement le fichier d'index .bai synchronisé avec le BAM dédupliqué.
        --VALIDATION_STRINGENCY SILENT
        # Ignore les avertissements bénins de format SAM sans interrompre le calcul.
    """
}
