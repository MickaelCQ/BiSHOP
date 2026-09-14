/*
    Reference: BiSHOP V2 CNV Consensus Engine
    Rôle : Fusion par chevauchement >= 50% des 3 callers CNV (CNVkit + ClinCNV + DECoN) avec tags CALLERS & NUM_CALLERS.
*/
process CONSENSUS_CNV {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/consensus_cnv" }, mode: 'copy'

    input:
    tuple val(meta), path(cnvkit_vcf), path(clincnv_tsv), path(decon_vcf)

    output:
    tuple val(meta), path("*.cnv_consensus.vcf.gz"), path("*.cnv_consensus.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    python3 ${projectDir}/bin/merge_cnv_consensus.py \\
        --cnvkit ${cnvkit_vcf} \\
        --clincnv ${clincnv_tsv} \\
        --decon ${decon_vcf} \\
        --sample ${meta.id} \\
        --out ${meta.id}.cnv_consensus.vcf

    bgzip -f ${meta.id}.cnv_consensus.vcf
    tabix -p vcf ${meta.id}.cnv_consensus.vcf.gz
    """
}
