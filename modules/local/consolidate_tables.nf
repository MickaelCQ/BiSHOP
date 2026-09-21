/*
    Reference: BiSHOP V2 Master Run Consolidator
    Rôle : Génération des deux fichiers CSV uniques de synthèse (SNV et CNV) dans ToBeConcluded/
*/
process CONSOLIDATE_RUN_TABLES {
    publishDir path: { "${params.outdir}/ToBeConcluded" }, mode: 'copy', overwrite: true

    input:
    path snv_csvs
    path clincnv_files
    path decon_csv
    path cnvkit_files

    output:
    path "Cohort_Master_SNV_INDEL.csv", emit: snv_master, optional: true
    path "Cohort_Master_CNV.csv",       emit: cnv_master, optional: true

    script:
    """
    python3 ${projectDir}/bin/consolidate_run_tables.py ${params.outdir} .
    """
}
