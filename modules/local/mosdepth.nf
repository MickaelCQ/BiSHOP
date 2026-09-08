/*
    Reference: mosdepth (Pedersen & Quinlan, 2018, Bioinformatics, DOI: 10.1093/bioinformatics/btx699)
*/
process MOSDEPTH {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.summary.txt"), emit: summary
    tuple val(meta), path("*.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.region.dist.txt"), emit: region_dist
    tuple val(meta), path("*.thresholds.bed.gz"), emit: thresholds

    script:
    """
    mosdepth \
        --by ${meta.bed} \
        --thresholds 1,10,20,30,50,100 \
        --fast-mode \
        --threads ${task.cpus} \
        ${meta.id} \
        ${bam}
    """
}
