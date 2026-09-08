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
    tuple val(meta), path("*.summary.txt"),     emit: summary
    tuple val(meta), path("*.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.region.dist.txt"), emit: region_dist
    tuple val(meta), path("*.thresholds.bed.gz"), emit: thresholds

    script:
    """
    mosdepth \\
        --by ${meta.bed} \\
        # Calcule la couverture moyenne et médiane restreinte aux régions codantes du BED d'activité du patient.
        --thresholds 1,10,20,30,50,100 \\
        # Paliers cliniques standard :
        #   >= 10X : seuil minimal de détection SNV.
        #   >= 30X : gold standard ACMG pour valider un variant constitutionnel hétérozygote.
        #   >= 100X : haute confiance pour les mosaïques et les panels constitutionnels.
        --fast-mode \\
        # Désactive le calcul de couverture base par base pour aller 10x plus vite en ne lisant que les bornes de reads.
        --threads ${task.cpus} \\
        ${meta.id} \\
        ${bam}
    """
}

