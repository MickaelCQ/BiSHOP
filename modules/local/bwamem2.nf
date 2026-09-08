/*
    Reference: bwa-mem2 (Vasimuddin et al., 2019, IEEE IPDPS, DOI: 10.1109/IPDPS.2019.00041)
*/
process BWAMEM2_ALIGN {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/aligned_bams" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path fasta_and_indexes

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_bai

    script:
    // Sélection automatique du fichier FASTA principal (hg38.fa)
    def fasta = fasta_and_indexes.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
    """
    bwa-mem2 mem -t ${task.cpus} -R '@RG\\tID:${meta.id}\\tSM:${meta.patient}\\tPL:ILLUMINA' ${fasta} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    """
}

