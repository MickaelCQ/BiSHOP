/*
    Reference: CNVkit (Talevich et al., 2016, PLOS Comput Biol, DOI: 10.1371/journal.pcbi.1004873)
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
    # 1. Batch coverage et segmentation Haar
    cnvkit.py batch ${bam} -n -f ${fasta} --targets ${meta.bed} --output-dir . -p ${task.cpus} --segment-method haar

    # 2. Appel des nombres de copies entières (méthode par seuil)
    cns_file=\$(ls *.cns | head -n 1)
    cnvkit.py call "\$cns_file" -m threshold -o ${meta.id}.call.cns

    # 3. Exportation au format VCF et compression
    cnvkit.py export vcf ${meta.id}.call.cns -o ${meta.id}.cnvkit.vcf
    bgzip -f ${meta.id}.cnvkit.vcf
    """
}
