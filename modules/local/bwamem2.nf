/*
    Reference: bwa-mem2 (Vasimuddin et al., 2019, IEEE IPDPS, DOI: 10.1101/IPDPS.2019.00041)
    Rôle : Alignement des reads appariés sur le génome de référence GRCh38 avec accélération vectorielle AVX-512.
*/
process BWAMEM2_ALIGN {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/aligned_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path fasta_and_indexes

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam_bai

    script:
    def fasta = fasta_and_indexes.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    # Exécution de BWA-MEM2 avec Read Group complet (traçabilité clinique stricte)
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        # Utilise l'ensemble des cœurs CPU alloués par SLURM.
        -R '@RG\\tID:${meta.id}\\tSM:${meta.patient}\\tPL:ILLUMINA\\tLB:WES_TWIST\\tPU:NextSeq2000' \\
        # @RG : Read Group obligatoire pour GATK :
        #   ID : Identifiant unique de la librairie.
        #   SM : Identifiant patient (pour le génotypage et la fusion de cohortes).
        #   PL : Plateforme de séquençage (déclenche les modèles d'erreur d'Illumina).
        #   LB : Type de librairie (WES).
        ${fasta} \\
        ${reads[0]} ${reads[1]} | \\
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        # Alloue 2 Go de RAM par thread pour le tri par coordonnées génomiques sans saturer la machine.
        -o ${meta.id}.sorted.bam -

    # Indexation immédiate du BAM (.bai) pour permettre les accès aléatoires rapides des callers
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    """
}
