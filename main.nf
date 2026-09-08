#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
========================================================================================
    PIPELINE CLINIQUE EXOME CONSTITUTIONNEL MULTI-ACTIVITES (MARFAN, SLA, PGX, TSA...)
========================================================================================
*/

include { FASTQ_ALIGN_WF     } from './subworkflows/FASTQ_ALIGN'
include { CALL_VARIANTS_WF   } from './subworkflows/CALL_VARIANTS'
include { MOSDEPTH           } from './modules/local/mosdepth'
include { BCFTOOLS_CONSENSUS } from './modules/local/bcftools_consensus'
include { VEP                } from './modules/local/vep'
include { MULTIQC            } from './modules/local/multiqc'

workflow {

    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def bed_path = row.bed_file ? file(row.bed_file) : file(params.default_bed)
            def meta = [
                id:       row.sample,
                patient:  row.patient,
                activity: row.activity,
                bed:      bed_path
            ]
            def fastq1 = row.fastq_1 ? file(row.fastq_1) : null
            def fastq2 = row.fastq_2 ? file(row.fastq_2) : null
            def bam    = row.bam ? file(row.bam) : null
            def bai    = row.bai ? file(row.bai) : null

            return [ meta, fastq1, fastq2, bam, bai ]
        }
        .branch { item ->
            def meta   = item[0]
            def fastq1 = item[1]
            def fastq2 = item[2]
            def bam    = item[3]
            def bai    = item[4]

            fastq: fastq1 != null
                return [ meta, [fastq1, fastq2] ]
            bam:   bam != null
                return [ meta, bam, bai ]
        }
        .set { ch_inputs }

    // 1. Alignement FASTQ + Déduplication Picard MarkDuplicates
    FASTQ_ALIGN_WF ( ch_inputs.fastq )

    // 2. Fusion des BAMs (ceux déjà prêts en entrée + ceux alignés et dédupliqués)
    ch_ready_bams = ch_inputs.bam.mix( FASTQ_ALIGN_WF.out.bam_bai )

    // 3. Mesure de couverture réelle (Mosdepth)
    MOSDEPTH ( ch_ready_bams )

    // 4. Variant Calling (DeepVariant, GATK HaplotypeCaller, FreeBayes, Manta, CNVkit)
    CALL_VARIANTS_WF ( ch_ready_bams )

    // 5. Consensus des 3 callers & Découpage BED Activité
    ch_snv_to_merge = CALL_VARIANTS_WF.out.dv_vcf
        .join( CALL_VARIANTS_WF.out.gatk_vcf, by: 0 )
        .join( CALL_VARIANTS_WF.out.freebayes_vcf, by: 0 )

    BCFTOOLS_CONSENSUS ( ch_snv_to_merge )

    // 6. Annotation VEP
    VEP ( BCFTOOLS_CONSENSUS.out.vcf_tbi, file(params.vep_cache), file(params.fasta) )

    // 7. Collecte et Rapport Global MultiQC (FastQC, Fastp, Picard Duplicates, Mosdepth, VEP)
    ch_qc_reports = FASTQ_ALIGN_WF.out.falco_txt.map{ meta, file -> file }
        .mix( FASTQ_ALIGN_WF.out.fastp_json.map{ meta, file -> file } )
        .mix( FASTQ_ALIGN_WF.out.markdup_metrics.map{ meta, file -> file } )
        .mix( MOSDEPTH.out.summary.map{ meta, file -> file } )
        .mix( MOSDEPTH.out.global_dist.map{ meta, file -> file } )
        .mix( VEP.out.summary.map{ meta, file -> file } )
        .collect()

    MULTIQC ( ch_qc_reports )
}
