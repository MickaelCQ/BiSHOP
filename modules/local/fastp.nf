/*
    Reference: fastp (Chen et al., 2018, Bioinformatics, DOI: 10.1093/bioinformatics/bty560)
    Rôle : Élimination adaptateurs, correction d'erreurs et rognage poly-G (NextSeq 2000).
*/
process FASTP {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/fastp/${meta.id}" }, mode: 'copy', pattern: '*.{html,json}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")         , emit: json
    tuple val(meta), path("*.html")         , emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_R1.trim.fastq.gz --out2 ${meta.id}_R2.trim.fastq.gz \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --qualified_quality_phred 20 \\
        --unqualified_percent_limit 30 \\
        --length_required 35 \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html
    """
}
