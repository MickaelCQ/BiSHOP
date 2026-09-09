/*
    Reference: FreeBayes (Garrison & Marth, 2012, arXiv:1207.3907)
    Rôle : Appel bayésien basé sur les micro-haplotypes.
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
    # Options FreeBayes :
    # -m 20 : Mapping Quality minimale >= 20
    # -q 20 : Base Quality minimale >= 20 (Q20)
    # --min-alternate-fraction 0.15 : VAF minimale >= 15% pour constitutionnel et mosaïques
    # --min-alternate-count 3 : Au moins 3 reads alternatifs
    # --min-coverage 10 : Profondeur minimale >= 10X
    freebayes \\
        -f ${fasta} \\
        -t ${meta.bed} \\
        -m 20 \\
        -q 20 \\
        --min-alternate-fraction 0.15 \\
        --min-alternate-count 3 \\
        --min-coverage 10 \\
        ${bam} | \\
    bgzip -c > ${meta.id}.freebayes.vcf.gz

    tabix -p vcf ${meta.id}.freebayes.vcf.gz
    """
}
