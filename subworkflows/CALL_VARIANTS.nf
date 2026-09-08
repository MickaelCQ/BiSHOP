include { DEEPVARIANT          } from '../modules/local/deepvariant'
include { GATK_HAPLOTYPECALLER } from '../modules/local/gatk_haplotypecaller'
include { FREEBAYES            } from '../modules/local/freebayes'
include { MANTA                } from '../modules/local/manta'
include { CNVKIT               } from '../modules/local/cnvkit'
include { CLINCNV              } from '../modules/local/clincnv'
include { DECON                } from '../modules/local/decon'

workflow CALL_VARIANTS_WF {
    take:
    ch_bams

    main:
    ch_fasta     = file(params.fasta)
    ch_fasta_fai = file(params.fasta_fai)
    ch_dict      = file(params.dict)

    // 1. SNV / Indels (Solo par échantillon)
    DEEPVARIANT ( ch_bams, ch_fasta, ch_fasta_fai )
    GATK_HAPLOTYPECALLER ( ch_bams, ch_fasta, ch_fasta_fai, ch_dict )
    FREEBAYES ( ch_bams, ch_fasta, ch_fasta_fai )

    // 2. Structural Variants (Solo par échantillon)
    MANTA ( ch_bams, ch_fasta, ch_fasta_fai )

    // 3. Copy Number Variants Solo (CNVkit avec référence plate)
    CNVKIT ( ch_bams, ch_fasta )

    // 4. Copy Number Variants Cohorte (DECoN & ClinCNV sur la série complète)
    ch_all_bams = ch_bams.map{ meta, bam, bai -> bam }.collect()
    ch_all_bais = ch_bams.map{ meta, bam, bai -> bai }.collect()
    ch_bed      = ch_bams.map{ meta, bam, bai -> meta.bed }.first()

    CLINCNV ( ch_all_bams, ch_all_bais, ch_bed, ch_fasta, ch_fasta_fai )
    DECON   ( ch_all_bams, ch_all_bais, ch_bed, ch_fasta, ch_fasta_fai )

    emit:
    dv_vcf        = DEEPVARIANT.out.vcf_tbi
    gatk_vcf      = GATK_HAPLOTYPECALLER.out.vcf_tbi
    freebayes_vcf = FREEBAYES.out.vcf_tbi
    manta_vcf     = MANTA.out.vcf_tbi
    cnv_vcf       = CNVKIT.out.vcf
    clincnv_tsv   = CLINCNV.out.tsv
    decon_csv     = DECON.out.csv
}
