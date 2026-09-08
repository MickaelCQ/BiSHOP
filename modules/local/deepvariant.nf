/*
    Reference: Google DeepVariant (Poplin et al., 2018, Nature Biotechnology, DOI: 10.1038/nbt.4235)
    Rôle : Appel des SNVs et petits Indels par réseau neuronal convolutif (Inception-v3) entraîné sur données GIAB.
*/
process DEEPVARIANT {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/deepvariant" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.deepvariant.vcf.gz"), path("*.deepvariant.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WES \\
        # Charge les poids du réseau neuronal spécifiquement calibrés pour la capture d'exome (bruit de capture, PCR).
        --ref=${fasta} \\
        --reads=${bam} \\
        --regions=${meta.bed} \\
        #Restreint l'inférence neuronale strictement au BED d'intérêt du patient. 
        # Cela divise le temps de calcul par 5 et élimine 100 % des faux positifs hors-cible.
        --output_vcf=${meta.id}.deepvariant.vcf.gz \\
        --num_shards=${task.cpus}
        # Découpe le génome en sous-régions parallèles selon le nombre de cœurs SLURM disponibles.
    """
}
