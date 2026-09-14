#!/usr/bin/env python3
"""
Compilateur de la Baseline Historique des Panels Custom (IonTorrent)
CHU de Nîmes — Norme ISO 15189
"""
import os
import glob
import pandas as pd

def compile_panel_baseline(annot_bed, pattern, out_csv, panel_name):
    print("=" * 80)
    print(f"COMPILATION BASELINE HISTORIQUE — PANEL : {panel_name}")
    print("=" * 80)
    
    # 1. Charger le BED annoté
    bed_df = pd.read_csv(annot_bed, sep="\t", header=None)
    bed_df = bed_df.iloc[:, [0, 1, 2, 3, 4, 5]]
    bed_df.columns = ["chr", "start", "end", "region_id", "gene", "type"]
    n_regions = len(bed_df)
    
    # 2. Chercher tous les fichiers Mosdepth correspondants
    qc_files = sorted(glob.glob(os.path.expanduser(f"~/Exome/BiSHOP_V2/legacy_panels/qc/*{pattern}*.regions.bed.gz")))
    print(f"Échantillons historiques trouvés : {len(qc_files)}")
    
    cov_cols = []
    for f in qc_files:
        sname = os.path.basename(f).replace("_legacy.regions.bed.gz", "")
        # Extraire la colonne de couverture
        cov_series = pd.read_csv(f, sep="\t", compression="gzip", header=None).iloc[:, -1]
        col_name = f"cov_{sname}"
        bed_df[col_name] = cov_series.values[:n_regions].round(1)
        cov_cols.append(col_name)
        print(f"  -> {sname} : profondeur moyenne = {bed_df[col_name].mean():.1f} X")
        
    if not cov_cols:
        print("Aucun fichier trouvé.")
        return
        
    # 3. Calculer les statistiques historiques par exon
    bed_df["Moyenne_Panel_Historique"] = bed_df[cov_cols].mean(axis=1).round(1)
    bed_df["Min_Panel_Historique"] = bed_df[cov_cols].min(axis=1).round(1)
    bed_df["Max_Panel_Historique"] = bed_df[cov_cols].max(axis=1).round(1)
    bed_df["Pct_Patients_ge_30X"] = (bed_df[cov_cols].apply(lambda x: (x >= 30).sum(), axis=1) / len(cov_cols) * 100).round(1)
    
    # Tri par gène et identifiant de région
    final_cols = ["gene", "region_id", "chr", "start", "end", "type", 
                  "Moyenne_Panel_Historique", "Min_Panel_Historique", "Max_Panel_Historique", "Pct_Patients_ge_30X"] + cov_cols
    final_df = bed_df[final_cols].sort_values(by=["gene", "region_id"])
    
    final_df.to_csv(out_csv, sep=";", index=False)
    
    # Bilan à l'écran
    glob_mean = bed_df["Moyenne_Panel_Historique"].mean()
    sous_30x = len(bed_df[bed_df["Moyenne_Panel_Historique"] < 30])
    
    print("-" * 80)
    print(f"Régions totales analysées                 : {n_regions}")
    print(f"Profondeur moyenne globale du Panel       : {glob_mean:.1f} X")
    print(f"Régions historiquement sous-couvertes (<30X) : {sous_30x} ({sous_30x/n_regions*100:.1f} %)")
    print(f"Fichier CSV de baseline généré             : {out_csv}")
    print("=" * 80 + "\n")

# Compilation pour SLA
compile_panel_baseline(
    "legacy_panels/bed/SLA_Panel_hg38_annotated.bed",
    "SLA",
    "legacy_panels/qc/SLA_legacy_baseline_coverage.csv",
    "SLA"
)

# Compilation pour TSA-HF
compile_panel_baseline(
    "legacy_panels/bed/TSA_Panel_hg38_annotated.bed",
    "TSA",
    "legacy_panels/qc/TSA_legacy_baseline_coverage.csv",
    "TSA-HF"
)
