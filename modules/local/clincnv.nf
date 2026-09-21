/*
    Reference: ClinCNV (Fugmann et al., 2020, Bioinformatics / Nature Communications, DOI: 10.1101/634048)
    Mode: Analyse multi-échantillons / Cohorte
*/
process CLINCNV {
    tag "Cohort (${bams.size()} samples)"
    publishDir path: { "${params.outdir}/vcfs/clincnv" }, mode: 'copy', overwrite: true

    input:
    path bams
    path bais
    path bed
    path fasta
    path fasta_fai

    output:
    path "*_cnvs.tsv",                   emit: tsv
    path "*_cnvs.seg",                   emit: seg, optional: true
    path "*.png",                        emit: png, optional: true
    path "ontargetNormal.summary.xls",   emit: xls, optional: true

    script:
    """
    export HOME=/tmp
    export R_LIBS_USER=""
    export PATH="/usr/local/bin:\$PATH"

    # 1. Construction du BED 5 colonnes (chr, start, end, GC, genes)
    /usr/local/bin/bedtools nuc -fi ${fasta} -bed ${bed} | \
        awk 'NR>1 {print \$1"\t"\$2"\t"\$3"\t"\$5"\t"\$4}' > clincnv_targets.bed

    # 2. Construction robuste de la ligne d'en-tête (100% Tabulations via printf)
    printf "chr\tstart\tend" > header.tsv
    for bam in ${bams}; do
        sname=\$(basename "\$bam" .markdup.bam | sed -e 's/\\.sorted//g' -e 's/\\.bam//g')
        printf "\t%s" "\$sname" >> header.tsv
    done
    printf "\n" >> header.tsv

    # 3. Construction des colonnes de données avec paste (tabulations)
    awk '{print \$1"\t"\$2"\t"\$3}' clincnv_targets.bed > cohort_cov.tmp
    for bam in ${bams}; do
        /usr/local/bin/bedtools coverage -a ${bed} -b "\$bam" -mean | \
            awk '{print \$NF}' > col.tmp
        paste cohort_cov.tmp col.tmp > cohort_cov.tmp2
        mv cohort_cov.tmp2 cohort_cov.tmp
    done

    # Assemblage final de la matrice
    cat header.tsv cohort_cov.tmp > cohort_normal.cov
    rm -f header.tsv cohort_cov.tmp col.tmp

    # 4. Exécution de ClinCNV
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
        --scoreG 20 \\
        --lengthG 2 \\
        --numberOfThreads ${task.cpus}

    # 5. Publication de TOUS les TSV, SEG, PNG et XLS
    find \${OUT_PATH}/normal/ -name "*_cnvs.tsv" -exec cp {} . \\;
    find \${OUT_PATH}/normal/ -name "*_cnvs.seg" -exec cp {} . \\;
    cp \${OUT_PATH}/ontargetNormal.summary.xls . 2>/dev/null || true
    cp \${OUT_PATH}/*.png . 2>/dev/null || true
    """
}
