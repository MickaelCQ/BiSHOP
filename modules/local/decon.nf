/*
    Reference: DECoN / ExomeDepth (Fowler et al., 2016, Wellcome Open Res, DOI: 10.12688/wellcomeopenres.10069.1)
    Mode: Analyse multi-échantillons / Cohorte (Requiert >= 3 échantillons)
*/
process DECON {
    tag "Cohort (${bams.size()} samples)"
    publishDir path: { "${params.outdir}/vcfs/decon" }, mode: 'copy'

    input:
    path bams
    path bais
    path bed
    path fasta
    path fasta_fai

    output:
    path "*.decon.csv", emit: csv, optional: true
    path "*.vcf.gz",    emit: vcf, optional: true

    script:
    def bams_list = bams.join(',')
    def bais_list = bais.join(',')
    """
    export HOME=/tmp
    export R_LIBS_USER=""

    # DECoN / ExomeDepth nécessite au moins 3 échantillons pour construire son aggregate reference set
    nb_samples=\$(echo "${bams_list}" | tr ',' '\n' | wc -l)
    if [ "\${nb_samples}" -lt 3 ]; then
        echo "================================================================================"
        echo "[INFO CLINIQUE] DECoN nécessite au moins 3 échantillons pour son modèle ExomeDepth."
        echo "Série actuelle : \${nb_samples} échantillons. Étape ignorée sans erreur."
        echo "================================================================================"
        echo -e "Sample,Chromosome,Start,End,Type,BF" > cohort.decon.csv
        exit 0
    fi

    # 1. Calcul du GC content par région BED avec bedtools
    echo "GC_CONTENT" > gc_content.tsv
    bedtools nuc -fi ${fasta} -bed ${bed} | awk 'NR>1 {print \$5}' >> gc_content.tsv

    # 2. Exécution de DECoN sur l'ensemble des BAMs de la cohorte
    Rscript /Rscript/DECoN.R ${bams_list} ${bed} ${bais_list} 0.01 gc_content.tsv
    mv DECoN_output.csv cohort.decon.csv 2>/dev/null || true

    # 3. Conversion de chaque échantillon au format VCF
    if [ -f cohort.decon.csv ] && [ -f sample_names.RData ]; then
        Rscript /Rscript/csv2vcf_DECoN.R cohort.decon.csv ${fasta} sample_names.RData
        for f in *.txt; do
            if [ -s "\$f" ]; then
                sample=\$(basename "\$f" .txt)
                cat "\$f" | bgzip -c > "\${sample}.decon.vcf.gz"
            fi
        done
    fi
    """
}
