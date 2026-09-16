/*
========================================================================================
    Module : Ensembl VEP 110.1 (Variant Effect Predictor)
    Référence : McLaren et al., 2016, Genome Biology, DOI: 10.1186/s13059-016-0974-4
    Rôle : Annotation fonctionnelle, clinique et structurale complète pour le diagnostic WES.
========================================================================================
*/
process VEP {
    tag "$meta.id ($meta.activity)"
    publishDir path: { "${params.outdir}/annotated_reports" }, mode: 'link', overwrite: true

    input:
    tuple val(meta), path(vcf), path(tbi)
    path vep_cache
    path fasta

    output:
    tuple val(meta), path("*.vep.vcf.gz"),  emit: vcf
    tuple val(meta), path("*_summary.txt"), emit: summary

    script:
    """
    # ==================================================================================
    # EXPLICATION CLINIQUE ET TECHNIQUE DES OPTIONS VEP DE BISHOP V2 :
    #
    # 1. ENTRÉES / SORTIES & PARALLÉLISATION :
    #    -i / -o           : Fichier VCF d'entrée et VCF de sortie compressé bgzip avec le
    #                        secteur d'activité officiel du patient (SLA, TSA_HF, MARFAN, PGX).
    #    --format vcf --vcf: Force l'entrée et la sortie en format VCF standardisé (RFC/GA4GH)
    #                        avec injection de toutes les annotations dans le champ INFO/CSQ.
    #    --fork 8          : Parallélisation multi-threadée sur les cœurs alloués par SLURM.
    #
    # 2. CACHE LOCAL & MODE HORS-LIGNE (OBLIGATOIRE CHU SANS INTERNET) :
    #    --dir_cache       : Emplacement du cache officiel Ensembl décompressé (homo_sapiens/110_GRCh38).
    #    --fasta           : Génome de référence GRCh38 pour vérifier les allèles et séquences protéiques.
    #    --offline         : Mode 100% hors-ligne interdisant toute requête HTTP vers l'extérieur.
    #    --assembly GRCh38 : Fixe l'assemblage génomique sur la version hg38.
    #
    # 3. LE SUPER-DRAPEAU CLINIQUE --everything (ACTIVATION DES STANDARDS INTERNATIONAUX) :
    #    --symbol          : Ajoute le symbole officiel du gène HUGO/HGNC (ex: FBN1, SOD1, FUS).
    #    --hgvs            : Calcule les nomenclatures internationales HGVS ADN (c.) et protéine (p.).
    #    --canonical       : Identifie le transcrit de référence médical de chaque gène.
    #    --mane            : Référence internationale MANE Select unifiée entre NCBI et Ensembl.
    #    --biotype         : Précise la nature biologique du transcrit (protein_coding, lncRNA, etc.).
    #    --numbers         : Numérote l'exon ou l'intron touché (ex: 3/15).
    #    --sift / --polyphen: Scores in silico historiques pour les variants faux-sens.
    #    --af_gnomad       : Fréquences alléliques globales gnomAD intégrées au cache (exomes et génomes).
    #
    # 4. BASES DE DONNÉES CLINIQUES EXTERNES (--custom & --plugin) :
    #    --custom ClinVar  : Fichier VCF NCBI à jour. Injecte :
    #                        - CLNSIG (Pathogenic, Likely_pathogenic, VUS, Benign)
    #                        - CLNREVSTAT (Niveau de preuve et nombre d'étoiles de confiance)
    #                        - CLNDN (Maladies et phénotypes associés MedGen)
    #    --custom AlphaMissense : Modèle IA DeepMind (2023) prédisant la déstabilisation 3D protéique :
    #                        - am_pathogenicity (score continu de 0 à 1)
    #                        - am_class (likely_pathogenic, ambiguous, likely_benign)
    #    --plugin REVEL    : Méta-score consensus ClinGen (Ioannidis et al.) combinant 13 prédicteurs
    #                        pour les variants faux-sens rares (seuil haute confiance > 0.75).
    #    --clin_sig_allele 1 : Associe la pathogénicité ClinVar strictement à l'allèle muté détecté.
    #
    # 5. STATISTIQUES POUR MULTIQC :
    #    --stats_file / --stats_text : Génère le rapport d'audit tabulé nécessaire à MultiQC.
    # ==================================================================================

    vep \\
        -i ${vcf} \\
        -o ${meta.id}.${meta.activity}.vep.vcf.gz \\
        --compress_output bgzip \\
        --format vcf --vcf \\
        --fork ${task.cpus} \\
        --dir_cache ${params.vep_cache} \\
        --fasta ${params.fasta} \\
        --offline \\
        --assembly GRCh38 \\
        --everything \\
        --custom ${params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
        --custom ${params.alphamissense},AlphaMissense,bed,exact,0,am_pathogenicity,am_class \\
        --plugin REVEL,file=${params.revel} \\
        --clin_sig_allele 1 \\
        --stats_file ${meta.id}.${meta.activity}_summary.txt \\
        --stats_text
    """
}
