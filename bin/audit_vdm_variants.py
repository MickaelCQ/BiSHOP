#!/usr/bin/env python3
"""
Audit Automatique de Validation de Méthode (VDM) — CHU de Nîmes
Scanne les VCFs des 7 callers sur les variants attendus
"""
import os
import subprocess
import pandas as pd

BASE_DIR = "results_nextseq_All_from_scratch/vcfs"

# Liste de vérité de vos variants cliniques connus
TRUTH_VARIANTS = [
    {
        "sample": "2603051686-120X",
        "gene": "PSEN1",
        "hgvs": "c.293_306delinsTG",
        "type": "SNV/INDEL",
        "region_hg38": "14:73170950-73171050"
    },
    {
        "sample": "2207131857",
        "gene": "SOD1",
        "hgvs": "c.358-10T>G",
        "type": "SNV/INDEL",
        "region_hg38": "21:31661350-31661450"
    },
    {
        "sample": "2107192220-A",
        "gene": "FUS",
        "hgvs": "c.1554_1557del",
        "type": "SNV/INDEL",
        "region_hg38": "16:31184300-31184400"
    },
    {
        "sample": "2006081825",
        "gene": "SHANK3",
        "hgvs": "c.3679dup",
        "type": "SNV/INDEL",
        "region_hg38": "22:50694700-50695100"
    },
    {
        "sample": "2408070565",
        "gene": "GRN",
        "hgvs": "Del exons 2 a 7",
        "type": "CNV",
        "region_hg38": "17:44345000-44355000"
    },
    {
        "sample": "2407011622",
        "gene": "LDLR",
        "hgvs": "Del exons 11 et 12",
        "type": "CNV",
        "region_hg38": "19:11110000-11125000"
    }
]

def query_vcf(vcf_path, region):
    if not os.path.exists(vcf_path):
        return False, ".", "."
    cmd = f"tabix {vcf_path} {region} 2>/dev/null"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = [l for l in res.stdout.strip().split("\n") if l]
    if lines:
        first = lines[0].split("\t")
        pos = first[1]
        ref = first[3]
        alt = first[4]
        qual = first[5]
        info = first[7]
        return True, qual, f"{ref}>{alt}"
    return False, ".", "."

def check_clincnv(sample, gene):
    tsv_path = f"{BASE_DIR}/clincnv/{sample}_cnvs.tsv"
    if not os.path.exists(tsv_path):
        return False, ".", "."
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        match = df[df["genes"].astype(str).str.contains(gene, na=False)]
        if len(match) > 0:
            cn = match.iloc[0]["CN_change"]
            qv = match.iloc[0]["qvalue"]
            return True, str(qv), f"CN={cn}"
    except Exception:
        pass
    return False, ".", "."

def check_decon(sample, gene):
    vcf_path = f"{BASE_DIR}/decon/{sample}.markdup.bam.decon.vcf.gz"
    if not os.path.exists(vcf_path):
        # Alternative naming
        vcf_path = f"{BASE_DIR}/decon/{sample}.decon.vcf.gz"
    if os.path.exists(vcf_path):
        cmd = f"zgrep '{gene}' {vcf_path} 2>/dev/null"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if res.stdout.strip():
            return True, "PASS", "DECoN_Call"
    return False, ".", "."

print("=" * 85)
print("AUDIT AUTOMATIQUE VDM : RECHERCHE DES MUTATIONS CLINIQUES ATTENDUES")
print("=" * 85)

results = []
for item in TRUTH_VARIANTS:
    s = item["sample"]
    g = item["gene"]
    reg = item["region_hg38"]
    
    # 1. SNV Callers
    dv_ok, dv_q, dv_call = query_vcf(f"{BASE_DIR}/deepvariant/{s}.deepvariant.vcf.gz", reg)
    gatk_ok, gatk_q, gatk_call = query_vcf(f"{BASE_DIR}/gatk/{s}.gatk.vcf.gz", reg)
    fb_ok, fb_q, fb_call = query_vcf(f"{BASE_DIR}/freebayes/{s}.freebayes.vcf.gz", reg)
    
    # 2. CNV Callers
    cnvkit_ok, cnvkit_q, cnvkit_call = query_vcf(f"{BASE_DIR}/cnvkit/{s}.cnvkit.vcf.gz", reg)
    clincnv_ok, clincnv_q, clincnv_call = check_clincnv(s, g)
    decon_ok, decon_q, decon_call = check_decon(s, g)
    
    results.append({
        "Patient": s,
        "Gene": g,
        "Variant_Attendu": item["hgvs"],
        "Type": item["type"],
        "DeepVariant": f"OUI ({dv_call})" if dv_ok else "NON",
        "GATK_HC": f"OUI ({gatk_call})" if gatk_ok else "NON",
        "FreeBayes": f"OUI ({fb_call})" if fb_ok else "NON",
        "ClinCNV": f"OUI ({clincnv_call})" if clincnv_ok else "NON",
        "DECoN": f"OUI ({decon_call})" if decon_ok else "NON",
        "CNVkit": "OUI" if cnvkit_ok else "NON"
    })

out_df = pd.DataFrame(results)
print(out_df[["Patient", "Gene", "Variant_Attendu", "DeepVariant", "GATK_HC", "FreeBayes", "ClinCNV", "DECoN"]].to_string(index=False))

out_df.to_csv("legacy_panels/vdm/Tableau_Audit_VDM_Auto.csv", sep=";", index=False)
print("\n" + "=" * 85)
print("Rapport VDM sauvegardé dans : legacy_panels/vdm/Tableau_Audit_VDM_Auto.csv")
print("=" * 85)
