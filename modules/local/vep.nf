/*
    Reference: Ensembl VEP (McLaren et al., 2016, Genome Biology, DOI: 10.1186/s13059-016-0974-4)
*/
process VEP {
    tag "$meta.id ($meta.activity)"
    publishDir path: { "${params.outdir}/annotated_reports" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path vep_cache
    path fasta

    output:
    tuple val(meta), path("*.vep.vcf.gz"),      emit: vcf
    tuple val(meta), path("*_summary.txt"),     emit: summary
    tuple val(meta), path("*_summary.html"),    emit: html, optional: true

    script:
    """
    vep \\
        -i ${vcf} \\
        -o ${meta.id}.${meta.activity}.vep.vcf.gz \\
        --compress_output bgzip \\
        --format vcf --vcf \\
        --fork ${task.cpus} \\
        --dir_cache ${vep_cache} \\
        --fasta ${fasta} \\
        --stats_file ${meta.id}.${meta.activity}_summary.txt \\
        --stats_text \\
        --hgvs --symbol --clin_sig_allele 1 --assembly GRCh38 --offline
    """
}
