include { FALCO          } from '../modules/local/falco'
include { FASTP          } from '../modules/local/fastp'
include { BWAMEM2_ALIGN  } from '../modules/local/bwamem2'
include { MARKDUPLICATES } from '../modules/local/markduplicates'

workflow FASTQ_ALIGN_WF {
    take:
    ch_fastq_inputs

    main:
    ch_fasta_and_indexes = files("${params.fasta}*")

    FALCO ( ch_fastq_inputs )
    FASTP ( ch_fastq_inputs )
    BWAMEM2_ALIGN ( FASTP.out.reads, ch_fasta_and_indexes )
    MARKDUPLICATES ( BWAMEM2_ALIGN.out.bam_bai )

    emit:
    bam_bai         = MARKDUPLICATES.out.bam_bai
    falco_txt       = FALCO.out.txt
    fastp_json      = FASTP.out.json
    markdup_metrics = MARKDUPLICATES.out.metrics
}
