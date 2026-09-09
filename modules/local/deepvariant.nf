/*
    Reference: Google DeepVariant (Poplin et al., 2018, Nature Biotechnology, DOI: 10.1038/nbt.4235)
    Rôle : Appel des SNVs et Indels par réseau neuronal convolutif entraîné sur données WES.
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
    # Options DeepVariant :
    # --model_type=WES : Modèle neuronal entraîné sur exomes (bruit de capture spécifique)
    # --regions : Restreint l'inférence strictement au BED d'intérêt (gain de temps x5)
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WES \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --regions=${meta.bed} \\
        --output_vcf=${meta.id}.deepvariant.vcf.gz \\
        --num_shards=${task.cpus}
    """
}
