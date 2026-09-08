/*
    Reference: CNVkit (Talevich et al., 2016, PLOS Comput Biol, DOI: 10.1371/journal.pcbi.1004873)
    Rôle : Détection de CNV (délétions/duplications d'exons) en mode échantillon unique sans cohorte.
*/
process CNVKIT {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/cnvkit" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.cnvkit.vcf.gz"), emit: vcf

    script:
    """
    # 1. Analyse par batch avec référence plate neutre
    cnvkit.py batch ${bam} \\
        -n \\
        # -n / --normal : Pas de panel de normaux requis, crée une référence virtuelle plate corrigée du GC.
        -f ${fasta} \\
        --targets ${meta.bed} \\
        --output-dir . \\
        -p ${task.cpus} \\
        --segment-method cbs
        # Utilise l'algorithme CBS (Circular Binary Segmentation), le standard le plus précis
        # en cytogénétique.

    # 2. Appel des nombres de copies discrets (Seuils log2 ratio constitutionnels)
    cns_file=\$(ls *.cns | head -n 1)
    cnvkit.py call "\$cns_file" \\
        -m threshold \\
        -t=-1.1,-0.4,0.3,0.7 \\
        # Seuils cliniques stricts de log2 ratio :
        #   < -1.1 : Délétion complète homozygote (CN=0)
        #   entre -1.1 et -0.4 : Délétion hétérozygote (CN=1)
        #   entre 0.3 et 0.7 : Duplication (CN=3)
        #   > 0.7 : Amplification majeure (CN>=4)
        --purity 1.0 \\
        # Échantillon constitutionnel pur à 100% (pas de contamination tumorale).
        -o ${meta.id}.call.cns

    # 3. Exportation au format VCF standardisé
    cnvkit.py export vcf ${meta.id}.call.cns -o ${meta.id}.cnvkit.vcf
    bgzip -f ${meta.id}.cnvkit.vcf
    """
}
