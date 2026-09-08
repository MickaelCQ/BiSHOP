/*
    Reference: GATK4 HaplotypeCaller (McKenna et al., 2010, Genome Res, DOI: 10.1101/gr.107524.110)
*/
process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/gatk" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.gatk.vcf.gz"), path("*.gatk.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    gatk HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}.gatk.vcf.gz \\
        -L ${meta.bed}
    """
}
