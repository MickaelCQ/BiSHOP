/*
    Reference: BCFtools (Li, 2011, Bioinformatics, DOI: 10.1093/bioinformatics/bt280)
    Rôle : Fusion des 3 callers SNV/Indels (DeepVariant + GATK + FreeBayes) et restriction au panel d'activité.
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
    # 1. Fusion des 3 callers SNV indépendants
    bcftools merge \\
        -m all \\
        # Fusionne les variants multi-alléliques sur une seule ligne VCF propre.
        --force-samples \\
        # Harmonise les identifiants d'échantillons identiques entre les différents VCFs.
        ${dv_vcf} ${gatk_vcf} ${fb_vcf} \\
        -O z -o merged.vcf.gz

    bcftools index -t merged.vcf.gz

    # 2. Découpage : Restriction aux ROI du panel d'activité (MARFAN, SLA, CARDIOMYOPATHIE, PGX...)
    bcftools view \\
        -R ${meta.bed} \\
        # -R : Filtre strictement les variants situés dans les exons du panel prescrit pour le patient.
        merged.vcf.gz \\
        -O z -o ${meta.id}.${meta.activity}.consolidated.activity.vcf.gz

    bcftools index -t ${meta.id}.${meta.activity}.consolidated.activity.vcf.gz
    """
}
