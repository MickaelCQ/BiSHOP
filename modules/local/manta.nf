/*
    Reference: Manta SV (Chen et al., 2016, Bioinformatics, DOI: 10.1093/bioinformatics/btw349)
    Rôle : Détection des délétions, duplications en tandem, inversions et translocations (> 50 pb).
*/
process MANTA {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/manta" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.manta.vcf.gz"), path("*.manta.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    # 1. Configuration de Manta
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        --exome \\
        # Désactive les modèles statistiques WGS basés sur la profondeur uniforme, 
        # et focalise Manta sur les reads discordants (chimeric reads) et split-reads aux bordures d'exons.
        --runDir manta_work

    # 2. Exécution parallèle
    manta_work/runWorkflow.py -j ${task.cpus}

    # 3. Extraction du VCF des variants constitutionnels diploïdes
    mv manta_work/results/variants/diploidSV.vcf.gz ${meta.id}.manta.vcf.gz
    mv manta_work/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.manta.vcf.gz.tbi
    """
}
