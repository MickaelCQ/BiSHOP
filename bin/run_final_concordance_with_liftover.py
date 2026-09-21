#!/usr/bin/env python3
import os, subprocess, pandas as pd
from pyliftover import LiftOver

BASE = os.path.expanduser("~/Exome/BiSHOP_V2")
WES_DIR = os.path.join(BASE, "result_11092026/aligned_bams")
LEGACY_DIR = os.path.join(BASE, "legacy_panels/bam")
OUT_DIR = os.path.join(BASE, "legacy_panels/concordance")
MOSDEPTH_DIR = os.path.join(OUT_DIR, "mosdepth_runs")
os.makedirs(MOSDEPTH_DIR, exist_ok=True)

# 1. RETRO-LIFTOVER VIA PYLIFTOVER (hg38 -> hg19)
print("=== 1. Traduction des coordonnées (hg38 -> hg19) ===")
chain_file = os.path.join(BASE, "legacy_panels/hg38ToHg19.over.chain.gz")
lo = LiftOver(chain_file)

def lift_bed_reverse(input_bed, output_bed):
    converted = 0
    with open(input_bed, 'r') as f_in, open(output_bed, 'w') as f_out:
        for line in f_in:
            p = line.strip().split('\t')
            if len(p) < 3: continue
            c, s, e, rest = p[0], int(p[1]), int(p[2]), p[3:]
            
            c_name = c if c.startswith('chr') else 'chr' + c
            ns = lo.convert_coordinate(c_name, s)
            ne = lo.convert_coordinate(c_name, e)
            
            if ns and ne and ns[0][0] == ne[0][0]:
                out_chr = ns[0][0].replace('chr', '')
                f_out.write(f"{out_chr}\t{ns[0][1]}\t{ne[0][1]}\t" + "\t".join(rest) + "\n")
                converted += 1
    print(f" -> {os.path.basename(input_bed)} : {converted} régions converties vers hg19.")

lift_bed_reverse(os.path.join(BASE, "legacy_panels/bed/SLA_Panel_GRCh38_curated.bed"), 
                 os.path.join(BASE, "legacy_panels/bed/SLA_Panel_hg19_curated.bed"))
lift_bed_reverse(os.path.join(BASE, "legacy_panels/bed/TSA_Panel_GRCh38_curated.bed"), 
                 os.path.join(BASE, "legacy_panels/bed/TSA_Panel_hg19_curated.bed"))

# 2. CONFIGURATION DES COHORTES (WES hg38 / Panels hg19)
COHORTS = {
    "SLA": {
        "bed_wes": "legacy_panels/bed/SLA_Panel_GRCh38_curated.bed",
        "bed_legacy": "legacy_panels/bed/SLA_Panel_hg19_curated.bed",
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
        "bed_legacy": "legacy_panels/bed/TSA_Panel_hg19_curated.bed",
        "pairs": [
            ("2006081825.markdup.bam", "2006081825_TSA-HF.bam", "2006081825"),
            ("2407011622.markdup.bam", "2407011622_TSA-HF.bam", "2407011622")
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

# 3. LANCEMENT MOSDEPTH ET FUSION
for panel, cfg in COHORTS.items():
    print(f"\n=== 2. Analyse Mosdepth & Fusion : {panel} ===")
    bed_wes, bed_legacy = os.path.join(BASE, cfg["bed_wes"]), os.path.join(BASE, cfg["bed_legacy"])
    
    # Le Squelette hg38 est le Maître
    df_master = pd.read_csv(bed_wes, sep='\t', header=None, names=['chr','start','end','region_id','gene','type'])
    wes_cols, panel_cols = [], []

    for wes_bam, leg_bam, sample_id in cfg["pairs"]:
        w_path, p_path = os.path.join(WES_DIR, wes_bam), os.path.join(LEGACY_DIR, leg_bam)
        w_pref, p_pref = os.path.join(MOSDEPTH_DIR, f"{panel}_WES_{sample_id}"), os.path.join(MOSDEPTH_DIR, f"{panel}_PANEL_{sample_id}")
        
        w_tsv = run_mosdepth(bed_wes, w_path, w_pref)
        p_tsv = run_mosdepth(bed_legacy, p_path, p_pref)
        
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
    df_master['Moyenne_WES_Cohorte'] = df_master[wes_cols].mean(axis=1).round(1) if wes_cols else 0.0
    df_master['Moyenne_Panel_Cohorte'] = df_master[panel_cols].mean(axis=1).round(1) if panel_cols else 0.0
    df_master['Ratio_Global_WES_Panel'] = df_master.apply(lambda r: round(r['Moyenne_WES_Cohorte'] / r['Moyenne_Panel_Cohorte'], 2) if r['Moyenne_Panel_Cohorte'] > 0 else 0.0, axis=1)
    df_master['Synthese_Accreditation'] = df_master['Moyenne_WES_Cohorte'].apply(get_status)
    
    final_cols = ['gene','region_id','chr','start','end','type','Synthese_Accreditation','Moyenne_WES_Cohorte','Moyenne_Panel_Cohorte','Ratio_Global_WES_Panel'] + wes_cols + panel_cols
    out_csv = os.path.join(OUT_DIR, f"{panel}_Master_Clinical_Concordance.csv")
    df_master[final_cols].to_csv(out_csv, sep=';', index=False)
    print(f"✅ Matrice validée et prête : {out_csv}")

