/*
    Reference: MultiQC (Ewels et al., 2016, Bioinformatics, DOI: 10.1093/bioinformatics/btw354)
*/
process MULTIQC {
    publishDir path: { "${params.outdir}/reports/multiqc" }, mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: html

    script:
    """
    multiqc . -f
    """
}
