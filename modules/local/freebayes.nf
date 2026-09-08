/*
    Reference: FreeBayes (Garrison & Marth, 2012, arXiv:1207.3907)
*/
process FREEBAYES {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/freebayes" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.freebayes.vcf.gz"), path("*.freebayes.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    freebayes \\
        -f ${fasta} \\
        -t ${meta.bed} \\
        -m 20 \\
        -q 20 \\
        ${bam} | bgzip -c > ${meta.id}.freebayes.vcf.gz

    tabix -p vcf ${meta.id}.freebayes.vcf.gz
    """
}
