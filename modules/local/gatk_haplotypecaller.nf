/*
    Reference: GATK4 HaplotypeCaller (McKenna et al., 2010, Genome Res, DOI: 10.1101/gr.107524.110)
    Rôle : Réassemblage local deNovo par graphe de de Bruijn et calcul de vraisemblance HMM.
*/
process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    publishDir path: { "${params.outdir}/vcfs/gatk" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.gatk.vcf.gz"), path("*.gatk.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga() - 4}g -XX:+UseParallelGC" HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}.gatk.vcf.gz \\
        -L ${meta.bed} \\
        # Cible uniquement les exons d'intérêt clinique du patient.
        -ip 100 \\
        # Interval Padding : Ajoute 100 paires de bases de part et d'autre de chaque exon.
        # Indispensable en génétique clinique pour capturer les variants d'épissage canoniques (+1/+2, -1/-2) 
        # et les variants introniques profonds régulateurs.
        --native-pair-hmm-threads ${task.cpus}
        # Parallélise le calcul mathématique le plus lourd de GATK (le modèle PairHMM en C++).
    """
}



