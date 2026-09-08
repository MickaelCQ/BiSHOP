/*
    Reference: fastp (Chen et al., 2018, Bioinformatics, DOI: 10.1093/bioinformatics/bty560)
    Rôle : Élimination des adaptateurs Illumina, correction d'erreurs et rognage des poly-G (NextSeq 2000).
*/
process FASTP {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/reports/fastp/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")         , emit: json
    tuple val(meta), path("*.html")         , emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_R1.trim.fastq.gz --out2 ${meta.id}_R2.trim.fastq.gz \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        # Détecte automatiquement les adaptateurs paired-end sans avoir à fournir de fichier FASTA d'adaptateurs.
        --trim_poly_g \\
        # Indispensable sur NextSeq 2000 (chimie 2 couleurs) : le 'G' correspond à l'absence de signal lumineux 
        # en fin de cycle (signal d'épuisement des réactifs) et doit être éliminé pour éviter les faux variants.
        --qualified_quality_phred 20 \\
        # Seuil de qualité Q20 (seules les bases ayant une précision >= 99.0% sont considérées de haute qualité).
        --unqualified_percent_limit 30 \\
        # Élimine les reads qui contiennent plus de 30% de bases de basse qualité (< Q20).
        --length_required 35 \\
        # Rejette les reads qui mesurent moins de 35 pb après rognage pour éviter les alignements non spécifiques.
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html
    """
}
