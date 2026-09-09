/*
    Reference: bwa-mem2 (Vasimuddin et al., 2019, IEEE IPDPS, DOI: 10.1101/IPDPS.2019.00041)
    Rôle : Alignement des reads appariés sur GRCh38 avec accélération vectorielle.
*/
process BWAMEM2_ALIGN {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/aligned_bams" }, mode: 'link'

    input:
    tuple val(meta), path(reads)
    path fasta_and_indexes

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam_bai

    script:
    def fasta = fasta_and_indexes.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    # Options BWA-MEM2 :
    # -R : Read Group complet obligatoire pour GATK (ID, SM=Patient, PL=Illumina, LB=WES, PU=NextSeq2000)
    # samtools sort -m 2G : Tri par coordonnées avec 2 Go de RAM par thread
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.patient}\\tPL:ILLUMINA\\tLB:WES_TWIST\\tPU:NextSeq2000' \\
        ${fasta} \\
        ${reads[0]} ${reads[1]} | \\
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o ${meta.id}.sorted.bam -

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    """
}
