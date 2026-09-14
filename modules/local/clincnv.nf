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

    # 2. Construction de la matrice de couverture de la cohorte
    awk '{print \$1"\t"\$2"\t"\$3}' clincnv_targets.bed > cohort_cov.tmp

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

    echo -e "\${sample_headers}" > cohort_normal.cov
    cat cohort_cov.tmp >> cohort_normal.cov
    rm -f cohort_cov.tmp *.col

    # 3. Exécution de ClinCNV
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

    # 4. Publication de TOUS les TSV, SEG, PNG et XLS
    find \${OUT_PATH}/normal/ -name "*_cnvs.tsv" -exec cp {} . \\;
    find \${OUT_PATH}/normal/ -name "*_cnvs.seg" -exec cp {} . \\;
    cp \${OUT_PATH}/ontargetNormal.summary.xls . 2>/dev/null || true
    cp \${OUT_PATH}/*.png . 2>/dev/null || true
    """
}
