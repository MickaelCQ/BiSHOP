/*
    Reference: Ensembl VEP (McLaren et al., 2016, Genome Biology, DOI: 10.1186/s13059-016-0974-4)
    Rôle : Annotation fonctionnelle, HGVS, transcrits canoniques et stats MultiQC.
*/
process VEP {
    tag "$meta.id ($meta.activity)"
    publishDir path: { "${params.outdir}/annotated_reports" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path vep_cache
    path fasta

    output:
    tuple val(meta), path("*.vep.vcf.gz"),  emit: vcf
    tuple val(meta), path("*_summary.txt"), emit: summary

    script:
    """
    # Options VEP :
    # --offline : Mode 100% hors-ligne sur cache local
    # --hgvs : Nomenclature HGVS officielle (c. et p.)
    # --symbol : Nom de gène HUGO
    # --canonical : Identifie le transcrit de référence
    # --clin_sig_allele 1 : Conserve la pathogénicité ClinVar
    # --stats_text : Rapport texte tabulé pour MultiQC
    vep \\
        -i ${vcf} \\
        -o ${meta.id}.${meta.activity}.vep.vcf.gz \\
        --compress_output bgzip \\
        --format vcf --vcf \\
        --fork ${task.cpus} \\
        --dir_cache ${vep_cache} \\
        --fasta ${fasta} \\
        --offline \\
        --assembly GRCh38 \\
        --hgvs \\
        --symbol \\
        --canonical \\
        --biotype \\
        --numbers \\
        --clin_sig_allele 1 \\
        --stats_file ${meta.id}.${meta.activity}_summary.txt \\
        --stats_text
    """
}
