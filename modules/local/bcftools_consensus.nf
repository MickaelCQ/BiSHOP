/*
    Reference: BiSHOP V2 Clinical Consensus Engine
    Rôle : Décomposition, alignement à gauche, consensus 3 callers avec tags CALLERS & NUM_CALLERS.
           Élimination native des sites de référence RefCall (NUM_CALLERS=0).
*/
process BCFTOOLS_CONSENSUS {
    tag "$meta.id ($meta.activity)"
    publishDir path: { "${params.outdir}/vcfs/consensus_activity" }, mode: 'copy'

    input:
    tuple val(meta), path(dv_vcf), path(dv_tbi), path(gatk_vcf), path(gatk_tbi), path(fb_vcf), path(fb_tbi)
    path fasta

    output:
    tuple val(meta), path("*.snv_consensus.vcf.gz"), path("*.snv_consensus.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    # 1. Normalisation, décomposition et élimination des non-variants (RefCall / AC=0)
    bcftools view -c 1 -e 'FILTER ~ "RefCall"' ${dv_vcf} | \
        bcftools norm -m -any -f ${fasta} -O z -o dv.norm.vcf.gz
    bcftools index -t dv.norm.vcf.gz

    bcftools view -c 1 ${gatk_vcf} | \
        bcftools norm -m -any -f ${fasta} -O z -o gatk.norm.vcf.gz
    bcftools index -t gatk.norm.vcf.gz

    bcftools view -c 1 ${fb_vcf} | \
        bcftools norm -m -any -f ${fasta} -O z -o fb.norm.vcf.gz
    bcftools index -t fb.norm.vcf.gz

    # 2. Fusion brute des 3 callers
    bcftools merge -m all --force-samples dv.norm.vcf.gz gatk.norm.vcf.gz fb.norm.vcf.gz -O v -o merged.raw.vcf

    # 3. Injection des tags cliniques et élimination stricte des 0-callers
    awk -F'\t' -v OFS='\t' -v smp="${meta.id}" '
    BEGIN { csq_added = 0 }
    /^##/ {
        if (!csq_added && /^##INFO/) {
            print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\\"Callers soutenant le variant (DeepVariant, GATK, FreeBayes)\\">"
            print "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\\"Nombre de callers en accord (1 a 3)\\">"
            csq_added = 1
        }
        print
        next
    }
    /^#CHROM/ {
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,smp
        next
    }
    {
        dv = (\$10 !~ /^\\.\\/\\./ && \$10 != "." && \$10 !~ /^0\\/0/) ? "DeepVariant" : ""
        gatk = (\$11 !~ /^\\.\\/\\./ && \$11 != "." && \$11 !~ /^0\\/0/) ? "GATK" : ""
        fb = (\$12 !~ /^\\.\\/\\./ && \$12 != "." && \$12 !~ /^0\\/0/) ? "FreeBayes" : ""
        
        callers = ""
        n = 0
        if (dv != "") { callers = dv; n++ }
        if (gatk != "") { callers = (callers != "" ? callers "," gatk : gatk); n++ }
        if (fb != "") { callers = (callers != "" ? callers "," fb : fb); n++ }
        
        # GARDE-FOU CLINIQUE ABSOLU : Si n == 0, on jette la ligne immédiatement !
        if (n == 0) next;
        
        best_gt = (\$10 !~ /^\\.\\/\\./ && \$10 != "." && \$10 !~ /^0\\/0/) ? \$10 : ((\$11 !~ /^\\.\\/\\./ && \$11 != "." && \$11 !~ /^0\\/0/) ? \$11 : \$12)
        flt = (n >= 2) ? "PASS" : "SINGLETON"
        
        \$7 = flt
        \$8 = \$8 ";CALLERS=" callers ";NUM_CALLERS=" n
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,best_gt
    }' merged.raw.vcf | bcftools view -O z -o merged.vcf.gz

    bcftools index -t merged.vcf.gz

    # 4. Restriction clinique au BED d activite du patient
    bcftools view -R ${meta.bed} merged.vcf.gz -O z -o ${meta.id}.${meta.activity}.snv_consensus.vcf.gz
    bcftools index -t ${meta.id}.${meta.activity}.snv_consensus.vcf.gz
    """
}
