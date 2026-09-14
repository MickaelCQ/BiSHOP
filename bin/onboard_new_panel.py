#!/usr/bin/env python3
"""
Intégrateur Universel de Nouveaux Panels Historiques (Clé en main)
Détection hg19/hg38 -> LiftOver auto -> Annotation Exon/Gène
BiSHOP V2 - CHU de Nîmes (Norme ISO 15189)
"""
import os
import sys
import argparse
import pandas as pd
from pyliftover import LiftOver

def parse_args():
    parser = argparse.ArgumentParser(description="Onboarding universel d'un nouveau panel clinique")
    parser.add_argument("--bed", required=True, help="Fichier BED historique du panel (ex: MARFAN_Panel.bed)")
    parser.add_argument("--name", required=True, help="Nom de l'activité / pathologie (ex: MARFAN)")
    parser.add_argument("--refseq", default=os.path.expanduser("~/Explorations/bioresources/RefSeq_hg38_CDS_neutral_canonical.bed"),
                        help="BED RefSeq CDS hg38 de référence")
    parser.add_argument("--chain", default=os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels/hg19ToHg38.over.chain.gz"),
                        help="Fichier de chaîne LiftOver")
    return parser.parse_args()

def is_hg19(df):
    # Teste les coordonnées du chromosome 1 (TARDBP ou coordonnées typiques)
    chr1_regions = df[df[0].astype(str).isin(['1', 'chr1'])]
    if len(chr1_regions) > 0:
        # Si une coordonnée est > 248956422 (taille max chr1 hg38), ou si position TARDBP hg19
        max_pos = chr1_regions[2].max()
        if max_pos > 248956422:
            return True
        # Heuristique sur région test
        tardbp_hg19 = chr1_regions[(chr1_regions[1] >= 11070000) & (chr1_regions[2] <= 11080000)]
        if len(tardbp_hg19) > 0:
            return True
    return False

def main():
    args = parse_args()
    print("=" * 80)
    print(f"INTÉGRATION DU NOUVEAU PANEL CLINIQUE : {args.name}")
    print("=" * 80)
    
    # 1. Lecture du BED brut
    raw_df = pd.read_csv(args.bed, sep="\t", header=None)
    print(f"Régions cibles brutes lues : {len(raw_df)}")
    
    # 2. Détection du génome
    check_hg19 = is_hg19(raw_df)
    out_dir = os.path.dirname(args.bed) or "."
    
    if check_hg19:
        print(">>> Génome détecté : GRCh37 / hg19. Lancement du LiftOver automatique vers hg38...")
        lo = LiftOver(args.chain)
        lifted_rows = []
        for _, row in raw_df.iterrows():
            c, s, e = str(row[0]), int(row[1]), int(row[2])
            c_name = c if c.startswith('chr') else 'chr' + c
            ns = lo.convert_coordinate(c_name, s)
            ne = lo.convert_coordinate(c_name, e)
            if ns and ne and ns[0][0] == ne[0][0]:
                out_c = ns[0][0].replace('chr', '')
                lifted_rows.append([out_c, ns[0][1], ne[0][1]])
        hg38_df = pd.DataFrame(lifted_rows)
        hg38_bed_path = os.path.join(out_dir, f"{args.name}_Panel_hg38.bed")
        hg38_df.to_csv(hg38_bed_path, sep="\t", header=False, index=False)
        print(f">>> LiftOver terminé : {len(hg38_df)}/{len(raw_df)} régions converties dans {hg38_bed_path}")
    else:
        print(">>> Génome détecté : GRCh38 / hg38. Pas de conversion requise.")
        hg38_df = raw_df[[0, 1, 2]]
        hg38_bed_path = args.bed

    # 3. Annotation automatique avec RefSeq CDS
    print(">>> Annotation des gènes et numérotation des exons...")
    ref_df = pd.read_csv(args.refseq, sep="\t", header=None)
    ref_df[0] = ref_df[0].astype(str).str.replace('chr', '')
    hg38_df[0] = hg38_df[0].astype(str).str.replace('chr', '')
    
    # Intersection simple en mémoire
    annotated_rows = []
    for _, row in hg38_df.iterrows():
        c, s, e = str(row[0]), int(row[1]), int(row[2])
        match = ref_df[(ref_df[0] == c) & (ref_df[1] < e) & (ref_df[2] > s)]
        if len(match) > 0:
            gene = match.iloc[0][3] if match.shape[1] >= 4 else f"Gene_{c}_{s}"
            annotated_rows.append([c, s, e, gene, "CDS"])
        else:
            annotated_rows.append([c, s, e, f"Target_{c}_{s}", "Non_CDS_Custom"])
            
    final_df = pd.DataFrame(annotated_rows, columns=["chr", "start", "end", "gene", "type"])
    final_df = final_df.sort_values(by=["chr", "start", "end"]).drop_duplicates(subset=["chr", "start", "end"])
    final_df["exon_num"] = final_df.groupby("gene").cumcount() + 1
    final_df["region_id"] = final_df["gene"] + "_Exon_" + final_df["exon_num"].astype(str) + "_" + final_df["type"]
    
    annot_bed_path = os.path.join(out_dir, f"{args.name}_Panel_hg38_annotated.bed")
    final_df[["chr", "start", "end", "region_id", "gene", "type"]].to_csv(annot_bed_path, sep="\t", header=False, index=False)
    
    print("=" * 80)
    print(f"PANEL {args.name} ENRICHI ET PRÊT POUR L'EXOME !")
    print(f"Fichier final : {annot_bed_path}")
    print("=" * 80 + "\n")

if __name__ == "__main__":
    main()
