#!/usr/bin/env python3
import os, subprocess, pandas as pd

BASE = os.path.expanduser("~/Exome/BiSHOP_V2")
WES_DIR = os.path.join(BASE, "result_11092026/aligned_bams")
LEGACY_DIR = os.path.join(BASE, "legacy_panels/bam")
OUT_DIR = os.path.join(BASE, "legacy_panels/concordance")
MOSDEPTH_DIR = os.path.join(OUT_DIR, "mosdepth_runs")
os.makedirs(MOSDEPTH_DIR, exist_ok=True)
REF_CDS = '/mnt/ngs_miseq/BIOINFO/Bioresources/beds/RefSeq_hg38_CDS_neutral_canonical.bed'

# 1. DICTIONNAIRE EXHAUSTIF DES GÈNES SUR LE BRIN NÉGATIF (GRCh38)
MINUS_STRAND_GENES = {
    'ANG', 'ANXA11', 'APP', 'CHCHD10', 'CHMP2B', 'FUS', 'HNRNPA1', 'OPTN', 'SQSTM1', 'TARDBP', # SLA
    'ABCB1', 'CHD8', 'CYP3A4', 'CYP3A5', 'DYRK1A', 'GRIN2B', 'MECP2', 'MED13L', 'NRXN1', # TSA
    'CYP2D6', 'MTHFR', 'NUDT15', 'UGT1A1', 'UGT1A4', 'UGT2B7' # PGX
}

# 2. COHORTE EXHAUSTIVE DE PRODUCTION
COHORTS = {
    "SLA": {
        "bed_raw": "legacy_panels/bed/SLA_Panel_hg38.bed",
        "bed_curated": "legacy_panels/bed/SLA_Panel_GRCh38_curated.bed",
        "pairs": [
            ("2107192220-A.markdup.bam", "2107192220_SLA.bam", "2107192220"),
            ("2207131857.markdup.bam", "2207131857_SLA.bam", "2207131857"),
            ("2408070565.markdup.bam", "2408070565_SLA.bam", "2408070565"),
            ("2409131379.markdup.bam", "2409131379_SLA.bam", "2409131379"),
            ("2603051686-120X.markdup.bam", "2603051686_SLA.bam", "2603051686"),
            ("ADN_CORIELL_NA24385-A.markdup.bam", "NA24385-2508121471_SLA.bam", "NA24385")
        ]
    },
    "TSA": {
        "bed_raw": "legacy_panels/bed/TSA_Panel_hg38.bed",
        "bed_curated": "legacy_panels/bed/TSA_Panel_GRCh38_curated.bed",
        "pairs": [
            ("2006081825.markdup.bam", "2006081825_TSA-HF.bam", "2006081825"),
            ("2407011622.markdup.bam", "2407011622_TSA-HF.bam", "2407011622")
        ]
    },
    "PGX": {
        "bed_raw": "legacy_panels/bed/PGX_Panel.bed",
        "bed_curated": "legacy_panels/bed/PGX_Panel_GRCh38_curated.bed",
        "pairs": [
            ("2312070839.markdup.bam", "2312070839_PGX.bam", "2312070839")
        ]
    }
}

def curate_bed(raw_bed, out_bed, panel_name):
    print(f"\n[1/3] Correction Biologique & Numérotation du Panel : {panel_name}")
    int_tmp, non_tmp = f"{out_bed}.intersect", f"{out_bed}.noncds"
    subprocess.run(f"bedtools intersect -wo -a {raw_bed} -b {REF_CDS} > {int_tmp}", shell=True)
    subprocess.run(f"bedtools intersect -v -a {raw_bed} -b {REF_CDS} > {non_tmp}", shell=True)

    df_cds = pd.read_csv(int_tmp, sep='\t', header=None)[[0, 1, 2, 6]].drop_duplicates()
    df_cds.columns = ['chr', 'start', 'end', 'gene']
    df_cds['type'] = 'CDS'

    df_non = pd.read_csv(non_tmp, sep='\t', header=None)[[0, 1, 2]]
    df_non.columns = ['chr', 'start', 'end']
    
    def resolve_non_cds(r):
        c, s = str(r['chr']).replace('chr', ''), int(r['start'])
        if c == '21' and 31659000 <= s <= 31669000: return 'SOD1', 'Intron_Custom'
        elif c == '1' and 11012000 <= s <= 11013000: return 'TARDBP', '5UTR_Regulateur'
        elif c == 'X' and 145800000 <= s <= 145810000: return 'AR', 'Regulateur'
        else: return f'Region_{c}_{s}', 'Non_CDS_Custom'

    assigned = df_non.apply(resolve_non_cds, axis=1)
    df_non['gene'] = [a[0] for a in assigned]
    df_non['type'] = [a[1] for a in assigned]

    full_df = pd.concat([df_cds, df_non], ignore_index=True)
    full_df['clean_chr'] = full_df['chr'].astype(str).str.replace('chr', '')
    
    curated_rows = []
    for gene, group in full_df.groupby('gene'):
        # Logique de tri ISO : Croissant si +, Décroissant si -
        is_plus = gene not in MINUS_STRAND_GENES
        group = group.sort_values(by=['start'], ascending=is_plus)
        
        cds_i, int_i = 1, 1
        for _, row in group.iterrows():
            rtype = row['type']
            if rtype == 'CDS':
                reg_id = f"{gene}_Exon_{cds_i}_CDS"
                cds_i += 1
            elif 'Intron' in rtype:
                reg_id = f"{gene}_Intron_{int_i}_Custom"
                int_i += 1
            else:
                reg_id = f"{gene}_{rtype}"
            curated_rows.append({'chr': row['clean_chr'], 'start': row['start'], 'end': row['end'], 'region_id': reg_id, 'gene': gene, 'type': rtype})

    df_curated = pd.DataFrame(curated_rows).sort_values(by=['chr', 'start'])
    df_curated[['chr', 'start', 'end', 'region_id', 'gene', 'type']].to_csv(out_bed, sep='\t', header=False, index=False)
    os.remove(int_tmp); os.remove(non_tmp)
    return df_curated

def run_mosdepth(bed, bam, prefix):
    if not os.path.exists(bam): return None
    subprocess.run(f"mosdepth -t 4 -b {bed} -n {prefix} {bam}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    tsv = f"{prefix}.tsv"
    subprocess.run(f"zcat {prefix}.regions.bed.gz > {tsv}", shell=True)
    return tsv

def get_status(c):
    if c >= 100: return "VALIDÉ_OPTIMAL (>=100X)"
    elif c >= 30: return "VALIDÉ_CONFORME (>=30X)"
    else: return "ALERTE_DROPOUT (<30X)"

for panel, cfg in COHORTS.items():
    raw_bed = os.path.join(BASE, cfg["bed_raw"])
    curated_bed = os.path.join(BASE, cfg["bed_curated"])
    
    if not os.path.exists(raw_bed): continue
    
    # Etape 1 : Nettoyage BED
    df_master = curate_bed(raw_bed, curated_bed, panel)
    wes_cols, panel_cols = [], []

    print(f"[2/3] Mosdepth sur la cohorte {panel} ({len(cfg['pairs'])} patients)")
    for wes_bam, leg_bam, sample_id in cfg["pairs"]:
        
        w_path, p_path = os.path.join(WES_DIR, wes_bam), os.path.join(LEGACY_DIR, leg_bam)
        w_pref, p_pref = os.path.join(MOSDEPTH_DIR, f"{panel}_WES_{sample_id}"), os.path.join(MOSDEPTH_DIR, f"{panel}_PANEL_{sample_id}")
        
        w_tsv = run_mosdepth(curated_bed, w_path, w_pref)
        p_tsv = run_mosdepth(curated_bed, p_path, p_pref)
        
        if w_tsv:
            df_master[f'cov_wes_{sample_id}'] = pd.read_csv(w_tsv, sep='\t', header=None)[4]
            wes_cols.append(f'cov_wes_{sample_id}')
        if p_tsv:
            df_master[f'cov_panel_{sample_id}'] = pd.read_csv(p_tsv, sep='\t', header=None)[4]
            panel_cols.append(f'cov_panel_{sample_id}')

    print(f"[3/3] Génération Matrice d'Accréditation : {panel}")
    df_master['Moyenne_WES_Cohorte'] = df_master[wes_cols].mean(axis=1).round(1) if wes_cols else 0.0
    df_master['Moyenne_Panel_Cohorte'] = df_master[panel_cols].mean(axis=1).round(1) if panel_cols else 0.0
    
    df_master['Ratio_Global_WES_Panel'] = df_master.apply(
        lambda r: round(r['Moyenne_WES_Cohorte'] / r['Moyenne_Panel_Cohorte'], 2) if r['Moyenne_Panel_Cohorte'] > 0 else 0.0, axis=1
    )
    df_master['Synthese_Accreditation'] = df_master['Moyenne_WES_Cohorte'].apply(get_status)
    
    final_cols = ['gene','region_id','chr','start','end','type','Synthese_Accreditation','Moyenne_WES_Cohorte','Moyenne_Panel_Cohorte','Ratio_Global_WES_Panel'] + wes_cols + panel_cols
    out_csv = os.path.join(OUT_DIR, f"{panel}_Master_Clinical_Concordance.csv")
    df_master[final_cols].to_csv(out_csv, sep=';', index=False)
    print(f"✅ Matrice Validée et Sauvegardée : {out_csv}\n")

