/*
    Reference: ClinCNV (Fugmann et al., 2020, Bioinformatics / Nature Communications, DOI: 10.1101/634048)
    Mode: Analyse multi-échantillons / Cohorte
*/
process CLINCNV {
    tag "Cohort (${bams.size()} samples)"
    publishDir path: { "${params.outdir}/vcfs/clincnv" }, mode: 'copy'

    input:
    path bams
    path bais
    path bed
    path fasta
    path fasta_fai

    output:
    path "*.clincnv.tsv", emit: tsv, optional: true
    path "*.vcf.gz",      emit: vcf, optional: true

    script:
    """
    export HOME=/tmp
    export R_LIBS_USER=""
    export PATH="/usr/local/bin:\$PATH"

    # 1. Préparation du BED strict à 5 colonnes (chr start end GC genes)
    /usr/local/bin/bedtools nuc -fi ${fasta} -bed ${bed} | \
        awk 'NR>1 {print \$1"\t"\$2"\t"\$3"\t"\$5"\t"\$4}' > clincnv_targets.bed

    # 2. Construction de la matrice de couverture multi-échantillons (Cohorte)
    # Initialisation avec les colonnes chr, start, end
    awk '{print \$1"\t"\$2"\t"\$3}' clincnv_targets.bed > cohort_cov.tmp

    # Ajout de la colonne de couverture pour chaque BAM de la cohorte
    sample_headers="chr\tstart\tend"
    for bam in ${bams}; do
        sname=\$(basename "\$bam" .markdup.bam)
        sname=\$(basename "\$sname" .sorted.bam)
        sname=\$(basename "\$sname" .bam)
        sample_headers="\${sample_headers}\t\${sname}"
        
        /usr/local/bin/bedtools coverage -a ${bed} -b "\$bam" -mean | \
            awk '{print \$NF}' > "\${sname}.col"
        
        paste cohort_cov.tmp "\${sname}.col" > cohort_cov.tmp2
        mv cohort_cov.tmp2 cohort_cov.tmp
    done

    # Ajout de l'en-tête officiel
    echo -e "\${sample_headers}" > cohort_normal.cov
    cat cohort_cov.tmp >> cohort_normal.cov
    rm -f cohort_cov.tmp *.col

    # 3. Exécution de ClinCNV sur la cohorte complète
    COV_PATH=\$(pwd)/cohort_normal.cov
    BED_PATH=\$(pwd)/clincnv_targets.bed
    OUT_PATH=\$(pwd)/clincnv_cohort_out

    Rscript /opt/clincnv/clinCNV.R \\
        --normal \${COV_PATH} \\
        --bed \${BED_PATH} \\
        --out \${OUT_PATH} \\
        --folderWithScript /opt/clincnv \\
        --colNum 4 \\
        --hg38 \\
        --numberOfThreads ${task.cpus}

    # 4. Récupération des résultats TSV et VCF de la cohorte
    if ls \${OUT_PATH}/*_cnvs.tsv 1> /dev/null 2>&1; then
        cp \${OUT_PATH}/*_cnvs.tsv .
    fi
    for f in \${OUT_PATH}/*.vcf; do
        if [ -s "\$f" ]; then
            bgzip -c "\$f" > "\$(basename "\$f").gz"
            tabix -p vcf "\$(basename "\$f").gz" 2>/dev/null || true
        fi
    done
    """
}
