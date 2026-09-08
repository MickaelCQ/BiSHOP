/*
    Reference: BCFtools (Li, 2011, Bioinformatics, DOI: 10.1093/bioinformatics/bt280)
*/
process BCFTOOLS_CONSENSUS {
    tag "$meta.id ($meta.activity)"
    publishDir path: { "${params.outdir}/vcfs/consensus_activity" }, mode: 'copy'

    input:
    tuple val(meta), path(dv_vcf), path(dv_tbi), path(gatk_vcf), path(gatk_tbi), path(fb_vcf), path(fb_tbi)

    output:
    tuple val(meta), path("*.consolidated.activity.vcf.gz"), path("*.consolidated.activity.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    # 1. Merge des 3 callers SNV
    bcftools merge -m all --force-samples ${dv_vcf} ${gatk_vcf} ${fb_vcf} -O z -o merged.vcf.gz
    bcftools index -t merged.vcf.gz

    # 2. Restriction aux ROI du BED de l'activité du patient (MARFAN, SLA, PGX...)
    bcftools view -R ${meta.bed} merged.vcf.gz -O z -o ${meta.id}.${meta.activity}.consolidated.activity.vcf.gz
    bcftools index -t ${meta.id}.${meta.activity}.consolidated.activity.vcf.gz
    """
}
