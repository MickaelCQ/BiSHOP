/*
    Reference: Manta SV (Chen et al., 2016, Bioinformatics, DOI: 10.1093/bioinformatics/btw349)
    Rôle : Détection des variants structuraux (SV > 50 pb).
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
    # Options Manta :
    # --exome : Focalise Manta sur les reads discordants et split-reads aux bordures d'exons
    # -j : Nombre de cœurs CPU alloués par SLURM
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        --exome \\
        --runDir manta_work

    manta_work/runWorkflow.py -j ${task.cpus}

    mv manta_work/results/variants/diploidSV.vcf.gz ${meta.id}.manta.vcf.gz
    mv manta_work/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.manta.vcf.gz.tbi
    """
}
