/*
    Reference: fastp (Chen et al., 2018, Bioinformatics, DOI: 10.1093/bioinformatics/bty560)
*/
process FASTP {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/fastp/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")         , emit: json

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_R1.trim.fastq.gz --out2 ${meta.id}_R2.trim.fastq.gz \\
        --trim_poly_g \\
        --thread ${task.cpus} \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html
    """
}
