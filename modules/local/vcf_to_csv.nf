process VCF_TO_CSV {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/clinical_csv" }, mode: 'link', overwrite: true

    input:
    tuple val(meta), path(vep_vcf)

    output:
    path "*.csv", emit: csv

    script:
    """
    python3 ${projectDir}/modules/local/vcf_explore_clinical.py \\
        --vcf ${vep_vcf} \\
        --sample ${meta.id} \\
        --activity ${meta.activity} \\
        --out ${meta.id}_${meta.activity}_clinical.csv
    """
}
