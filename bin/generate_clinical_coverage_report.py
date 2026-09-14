#!/usr/bin/env python3
"""
Rapport de Couverture Comparative Clinique par Gène et par Exon
Panel Custom (Legacy IonTorrent) vs Exome WES (BiSHOP V2 Illumina)
Laboratoire de Génétique Médicale — CHU de Nîmes (Norme ISO 15189)
"""
import sys
import gzip
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Générateur de rapport clinique comparatif par Gène/Exon")
    parser.add_argument("--legacy-regions", required=True, help="Sortie Mosdepth *.regions.bed.gz de l'ancien Panel")
    parser.add_argument("--wes-regions", required=True, help="Sortie Mosdepth *.regions.bed.gz du nouvel Exome WES")
    parser.add_argument("--annot-bed", required=True, help="BED annoté officiel hg38 (SLA_Panel_hg38_annotated.bed ou TSA)")
    parser.add_argument("--sample", required=True, help="Identifiant du patient (ex: 2207131857)")
    parser.add_argument("--min-cov", type=float, default=30.0, help="Seuil de couverture minimale diagnostique (defaut: 30X)")
    parser.add_argument("--out-csv", required=True, help="Chemin du fichier CSV de sortie")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Chargement du BED annoté officiel
    annot_df = pd.read_csv(args.annot_bed, sep="\t", header=None)
    annot_df = annot_df.iloc[:, [0, 1, 2, 3, 4, 5]]
    annot_df.columns = ["chr", "start", "end", "region_id", "gene", "type"]
    
    # 2. Chargement de Mosdepth Legacy (Panel)
    leg_df = pd.read_csv(args.legacy_regions, sep="\t", compression="gzip", header=None)
    leg_cov = leg_df.iloc[:, -1].values
    
    # 3. Chargement de Mosdepth WES (Exome)
    wes_df = pd.read_csv(args.wes_regions, sep="\t", compression="gzip", header=None)
    wes_cov = wes_df.iloc[:, -1].values
    
    # Construction du tableau comparatif direct ligne par ligne
    n_regions = len(annot_df)
    annot_df["cov_legacy"] = pd.Series(leg_cov[:n_regions]).round(1)
    annot_df["cov_wes"] = pd.Series(wes_cov[:n_regions]).round(1)
    
    # Calcul du ratio et du statut clinique
    annot_df["ratio_wes_vs_panel"] = (annot_df["cov_wes"] / annot_df["cov_legacy"].replace(0, 0.1)).round(2)
    
    def get_status(row):
        if row["cov_wes"] >= 100:
            return "OPTIMAL (>=100X)"
        elif row["cov_wes"] >= args.min_cov:
            return "CONFORME (>=30X)"
        elif row["cov_legacy"] >= args.min_cov and row["cov_wes"] < args.min_cov:
            return "DROPOUT_CRITIQUE (<30X)"
        else:
            return "SOUS_COUVERTURE (<30X)"
            
    annot_df["Statut_Clinique"] = annot_df.apply(get_status, axis=1)
    
    final_cols = [
        "gene", "region_id", "chr", "start", "end", "type",
        "cov_legacy", "cov_wes", "ratio_wes_vs_panel", "Statut_Clinique"
    ]
    final_df = annot_df[final_cols].sort_values(by=["gene", "region_id"])
    final_df.to_csv(args.out_csv, sep=";", index=False)
    
    # Bilan résumé
    total = len(final_df)
    optimal = len(final_df[final_df["Statut_Clinique"].str.startswith("OPTIMAL")])
    conforme = len(final_df[final_df["Statut_Clinique"].str.startswith("CONFORME")])
    dropouts = len(final_df[final_df["Statut_Clinique"].str.startswith("DROPOUT")])
    
    print("=" * 80)
    print(f"RAPPORT CLINIQUE DE CONCORDANCE DE COUVERTURE — PATIENT : {args.sample}")
    print("=" * 80)
    print(f"Total des régions cibles évaluées : {total}")
    print(f"Régions optimales (>= 100X en WES): {optimal} ({optimal/total*100:.1f} %)")
    print(f"Régions conformes (>= 30X en WES) : {conforme + optimal} ({(conforme+optimal)/total*100:.1f} %)")
    print(f"ALERTES DROPOUTS (< 30X en WES)   : {dropouts} ({dropouts/total*100:.1f} %)")
    print("=" * 80)
    print(f"Fichier CSV médical généré : {args.out_csv}\n")

if __name__ == "__main__":
    main()
