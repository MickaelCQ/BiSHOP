#!/usr/bin/env python3
import os, subprocess, pandas as pd

BASE = os.path.expanduser("~/Exome/BiSHOP_V2")
WES_DIR = os.path.join(BASE, "result_11092026/aligned_bams")
LEGACY_DIR = os.path.join(BASE, "legacy_panels/bam")
OUT_DIR = os.path.join(BASE, "legacy_panels/concordance")
MOSDEPTH_DIR = os.path.join(OUT_DIR, "mosdepth_runs")

COHORTS = {
    "SLA": {
        "bed_wes": "legacy_panels/bed/SLA_Panel_GRCh38_curated.bed",
        "bed_legacy": "legacy_panels/bed/SLA_Panel_hg19_curated.bed", # BED hg19 !
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
        "bed_wes": "legacy_panels/bed/TSA_Panel_GRCh38_curated.bed",
        "bed_legacy": "legacy_panels/bed/TSA_Panel_hg19_curated.bed", # BED hg19 !
        "pairs": [
            ("2006081825.markdup.bam", "2006081825_TSA-HF.bam", "2006081825"),
            ("2407011622.markdup.bam", "2407011622_TSA-HF.bam", "2407011622")
        ]
    },
    "PGX": {
        "bed_wes": "legacy_panels/bed/PGX_Panel_GRCh38_curated.bed",
        "bed_legacy": "legacy_panels/bed/PGX_Panel_GRCh38_curated.bed", # PGX reste hg38 !
        "pairs": [
            ("2312070839.markdup.bam", "2312070839_PGX.bam", "2312070839")
        ]
    }
}

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
    print(f"\n[Calcul] Analyse Mosdepth & Fusion hg38/hg19 : {panel}")
    bed_wes, bed_legacy = os.path.join(BASE, cfg["bed_wes"]), os.path.join(BASE, cfg["bed_legacy"])
    
    # 1. Squelette hg38 = Maître absolu
    df_master = pd.read_csv(bed_wes, sep='\t', header=None, names=['chr','start','end','region_id','gene','type'])
    wes_cols, panel_cols = [], []

    for wes_bam, leg_bam, sample_id in cfg["pairs"]:
        w_path, p_path = os.path.join(WES_DIR, wes_bam), os.path.join(LEGACY_DIR, leg_bam)
        w_pref, p_pref = os.path.join(MOSDEPTH_DIR, f"{panel}_WES_{sample_id}"), os.path.join(MOSDEPTH_DIR, f"{panel}_PANEL_{sample_id}")
        
        w_tsv = run_mosdepth(bed_wes, w_path, w_pref)
        p_tsv = run_mosdepth(bed_legacy, p_path, p_pref)
        
        # 2. Merge sur nom de région (region_id) -> Les chiffres Panel rejoignent l'espace hg38 !
        if w_tsv:
            df_w = pd.read_csv(w_tsv, sep='\t', header=None, names=['chr','start','end','region_id','cov'])
            df_master = df_master.merge(df_w[['region_id', 'cov']], on='region_id', how='left')
            df_master.rename(columns={'cov': f'cov_wes_{sample_id}'}, inplace=True)
            wes_cols.append(f'cov_wes_{sample_id}')
            
        if p_tsv:
            df_p = pd.read_csv(p_tsv, sep='\t', header=None, names=['chr','start','end','region_id','cov'])
            df_master = df_master.merge(df_p[['region_id', 'cov']], on='region_id', how='left')
            df_master.rename(columns={'cov': f'cov_panel_{sample_id}'}, inplace=True)
            panel_cols.append(f'cov_panel_{sample_id}')

    df_master.fillna(0.0, inplace=True)
    
    # 3. Statistiques finales
    df_master['Moyenne_WES_Cohorte'] = df_master[wes_cols].mean(axis=1).round(1) if wes_cols else 0.0
    df_master['Moyenne_Panel_Cohorte'] = df_master[panel_cols].mean(axis=1).round(1) if panel_cols else 0.0
    df_master['Ratio_Global_WES_Panel'] = df_master.apply(lambda r: round(r['Moyenne_WES_Cohorte'] / r['Moyenne_Panel_Cohorte'], 2) if r['Moyenne_Panel_Cohorte'] > 0 else 0.0, axis=1)
    df_master['Synthese_Accreditation'] = df_master['Moyenne_WES_Cohorte'].apply(get_status)
    
    final_cols = ['gene','region_id','chr','start','end','type','Synthese_Accreditation','Moyenne_WES_Cohorte','Moyenne_Panel_Cohorte','Ratio_Global_WES_Panel'] + wes_cols + panel_cols
    out_csv = os.path.join(OUT_DIR, f"{panel}_Master_Clinical_Concordance.csv")
    df_master[final_cols].to_csv(out_csv, sep=';', index=False)
    print(f"✅ Matrice validée : {out_csv}")

