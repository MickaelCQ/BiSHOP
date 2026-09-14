#!/usr/bin/env python3
"""
Enrichissement Clinique des Panels Cibles hg38 : Noms de Gènes, Numéros d'Exons et Statut
BiSHOP V2 - CHU de Nîmes (Norme ISO 15189)
"""
import os
import sys
import pandas as pd

BASE_DIR = os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels")

def build_annotated_panel(intersect_file, non_cds_file, output_bed, panel_name):
    # 1. Charger les régions qui intersectent avec RefSeq CDS
    df_cds = pd.read_csv(intersect_file, sep="\t", header=None)
    # Colonnes bedtools -wo : 0=chr, 1=start, 2=end (panel) | 6=gene (RefSeq)
    cds_regions = df_cds[[0, 1, 2, 6]].drop_duplicates()
    cds_regions.columns = ["chr", "start", "end", "gene"]
    cds_regions["type"] = "CDS"
    
    # 2. Charger les régions hors CDS et leur assigner leur identité clinique
    try:
        df_non_cds = pd.read_csv(non_cds_file, sep="\t", header=None)
        df_non_cds = df_non_cds[[0, 1, 2]].drop_duplicates()
        df_non_cds.columns = ["chr", "start", "end"]
        
        def assign_gene(row):
            c, s = str(row['chr']), int(row['start'])
            if c in ['1', 'chr1'] and 11000000 <= s <= 11050000:
                return "TARDBP", "5UTR_Regulateur"
            elif c in ['21', 'chr21'] and 31650000 <= s <= 31670000:
                return "SOD1", "Exon_Intron_Custom"
            elif c in ['X', 'chrX'] and 145000000 <= s <= 146000000:
                return "AR", "Regulateur"
            else:
                return f"Region_{c}_{s}", "Non_CDS_Custom"
                
        assigned = df_non_cds.apply(assign_gene, axis=1)
        df_non_cds["gene"] = [a[0] for a in assigned]
        df_non_cds["type"] = [a[1] for a in assigned]
        
        full_df = pd.concat([cds_regions, df_non_cds], ignore_index=True)
    except Exception as e:
        print(f"[{panel_name}] Note sur non-CDS : {e}")
        full_df = cds_regions

    # Nettoyage et tri
    full_df["clean_chr"] = full_df["chr"].astype(str).str.replace("chr", "")
    full_df["start"] = full_df["start"].astype(int)
    full_df["end"] = full_df["end"].astype(int)
    full_df = full_df.sort_values(by=["clean_chr", "start", "end"]).drop_duplicates(subset=["clean_chr", "start", "end"])
    
    # 3. Numérotation automatique des exons par gène
    full_df["exon_num"] = full_df.groupby("gene").cumcount() + 1
    full_df["region_id"] = full_df["gene"] + "_Exon_" + full_df["exon_num"].astype(str) + "_" + full_df["type"]
    
    # Sauvegarde du BED final annoté (chr, start, end, region_id, gene, type)
    final_bed = full_df[["clean_chr", "start", "end", "region_id", "gene", "type"]]
    final_bed.to_csv(output_bed, sep="\t", header=False, index=False)
    
    print(f"[{panel_name}] Panel annoté généré avec succès : {len(final_bed)} régions dans {os.path.basename(output_bed)}")

# Exécution avec chemins absolus
build_annotated_panel(
    os.path.join(BASE_DIR, "intersection/SLA_hg38_intersect_RefSeq.bed"),
    os.path.join(BASE_DIR, "intersection/SLA_non_CDS_regions.bed"),
    os.path.join(BASE_DIR, "bed/SLA_Panel_hg38_annotated.bed"),
    "SLA"
)

build_annotated_panel(
    os.path.join(BASE_DIR, "intersection/TSA_hg38_intersect_RefSeq.bed"),
    os.path.join(BASE_DIR, "intersection/TSA_non_CDS_regions.bed"),
    os.path.join(BASE_DIR, "bed/TSA_Panel_hg38_annotated.bed"),
    "TSA"
)
