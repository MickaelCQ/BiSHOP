/*
    Reference: FreeBayes (Garrison & Marth, 2012, arXiv:1207.3907)
    Rôle : Appel bayésien basé sur les haplotypes courts, très sensible pour les Indels complexes.
*/
process FREEBAYES {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/freebayes" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.freebayes.vcf.gz"), path("*.freebayes.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    freebayes \\
        -f ${fasta} \\
        -t ${meta.bed} \\
        # Restreint l'analyse aux coordonnées du BED.
        -m 20 \\
        # Mapping Quality minimale >= 20 (élimine les reads alignés à plusieurs endroits du génome).
        -q 20 \\
        # Base Quality minimale >= 20 (Q20 = 99% de certitude sur le nucléotide).
        --min-alternate-fraction 0.20 \\
        # Fréquence allélique minimale (VAF >= 15%) pour détecter les hétérozygotes constitutionnels (autour de 50%) 
        # tout en capturant d'éventuelles mosaïques germinales sans bruit de fond.
        --min-alternate-count 5 \\
        # Exige au moins 3 reads indépendants portant le variant pour confirmer l'appel.
        --min-coverage 10 \\
        # Profondeur totale minimale de 10X au locus.
        ${bam} | \\
    bgzip -c > ${meta.id}.freebayes.vcf.gz

    # Indexation standard tabix
    tabix -p vcf ${meta.id}.freebayes.vcf.gz
    """
}
