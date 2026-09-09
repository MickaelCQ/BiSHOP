/*
    Reference: mosdepth (Pedersen & Quinlan, 2018, Bioinformatics, DOI: 10.1093/bioinformatics/btx699)
    Rôle : Mesure ultra-rapide de la profondeur de séquençage sur les régions cibles du patient.
*/
process MOSDEPTH {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/mosdepth" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.summary.txt"),       emit: summary
    tuple val(meta), path("*.global.dist.txt"),   emit: global_dist
    tuple val(meta), path("*.region.dist.txt"),   emit: region_dist
    tuple val(meta), path("*.thresholds.bed.gz"), emit: thresholds

    script:
    """
    # Options mosdepth :
    # --by : Calcule la couverture restreinte aux régions codantes du BED d'activité
    # --thresholds : Paliers cliniques 1, 10, 20, 30, 50, 100X (ACMG standard)
    # --fast-mode : Ne lit que les bornes de reads (10x plus rapide)
    mosdepth \\
        --by ${meta.bed} \\
        --thresholds 1,10,20,30,50,100 \\
        --fast-mode \\
        --threads ${task.cpus} \\
        ${meta.id} \\
        ${bam}
    """
}
