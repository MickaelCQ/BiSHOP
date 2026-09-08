/*
    Reference: DECoN / ExomeDepth (Fowler et al., 2016, Wellcome Open Res, DOI: 10.12688/wellcomeopenres.10069.1)
    Rôle : Détection fine des CNVs exoniques par agrégation des échantillons les plus corrélés de la série.
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

    # Garde-fou méthodologique : ExomeDepth a besoin d'au moins 3 échantillons pour bâtir son modèle de régression
    nb_samples=\$(echo "${bams_list}" | tr ',' '\n' | wc -l)
    if [ "\${nb_samples}" -lt 3 ]; then
        echo "[INFO CLINIQUE] DECoN nécessite >= 3 échantillons. Série actuelle : \${nb_samples}. Étape ignorée."
        echo -e "Sample,Chromosome,Start,End,Type,BF" > cohort.decon.csv
        exit 0
    fi

    # 1. Calcul du GC content par région BED
    echo "GC_CONTENT" > gc_content.tsv
    bedtools nuc -fi ${fasta} -bed ${bed} | awk 'NR>1 {print \$5}' >> gc_content.tsv

    # 2. Exécution de DECoN sur l'ensemble de la série
    Rscript /Rscript/DECoN.R \\
        ${bams_list} \\
        ${bed} \\
        ${bais_list} \\
        0.01 \\
        # trans_prob = 0.01 : Probabilité de transition du modèle de Markov caché (HMM) entre état normal et altéré.
        gc_content.tsv

    mv DECoN_output.csv cohort.decon.csv 2>/dev/null || true

    # 3. Conversion automatique au format VCF par échantillon
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
