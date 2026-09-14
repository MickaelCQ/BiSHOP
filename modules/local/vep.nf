/*
    Reference: Ensembl VEP (McLaren et al., 2016, Genome Biology, DOI: 10.1186/s13059-016-0974-4)
    Rôle : Annotation clinique de haute précision (conserve CALLERS, NUM_CALLERS, HGVS, gènes).
*/
process VEP {
    tag "$meta.id ($meta.var_type)"
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
    vep \\
        -i ${vcf} \\
        -o ${meta.id}.${meta.var_type}.vep.vcf.gz \\
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
        --stats_file ${meta.id}.${meta.var_type}_summary.txt \\
        --stats_text
    """
}
