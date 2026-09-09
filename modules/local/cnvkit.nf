/*
    Reference: CNVkit (Talevich et al., 2016, PLOS Comput Biol, DOI: 10.1371/journal.pcbi.1004873)
    Rôle : Détection de CNV solo avec référence plate neutre.
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
    # Options CNVkit :
    # -n : Référence plate sans panel de normaux
    # --segment-method cbs : Segmentation binaire circulaire
    # -t=-1.1,-0.4,0.3,0.7 : Seuils cliniques de log2 ratio (homozygote del, hétérozygote del, normal, dup, amp)
    # --purity 1.0 : Échantillon constitutionnel pur
    cnvkit.py batch ${bam} \\
        -n \\
        -f ${fasta} \\
        --targets ${meta.bed} \\
        --output-dir . \\
        -p ${task.cpus} \\
        --segment-method cbs

    cns_file=\$(ls *.cns | head -n 1)
    cnvkit.py call "\$cns_file" \\
        -m threshold \\
        -t=-1.1,-0.4,0.3,0.7 \\
        --purity 1.0 \\
        -o ${meta.id}.call.cns

    cnvkit.py export vcf ${meta.id}.call.cns -o ${meta.id}.cnvkit.vcf
    bgzip -f ${meta.id}.cnvkit.vcf
    """
}
