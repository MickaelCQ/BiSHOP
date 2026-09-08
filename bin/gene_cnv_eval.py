#!/usr/bin/env python3
"""
Evaluateur de CNV a l'echelle du gene (recouvrement >= 50%)
BiSHOP V2 - CHU de Nimes
"""
import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Evaluation clinique de CNV a l'echelle du gene (Overlap >= 50%)")
    parser.add_argument("--truth", required=True, help="Fichier de verite (Colonnes: chr, start, end, gene, type)")
    parser.add_argument("--calls", required=True, help="Fichier d'appels CNV (Colonnes: chr, start, end, type)")
    parser.add_argument("--overlap", type=float, default=0.50, help="Seuil de recouvrement minimum du gene (defaut: 0.50)")
    parser.add_argument("--out", default="cnv_eval_report.txt", help="Fichier texte de sortie")
    return parser.parse_args()

def normalize_chr(c):
    c = str(c).strip()
    return c if c.startswith("chr") else "chr" + c

def normalize_type(t):
    t = str(t).upper().strip()
    if any(x in t for x in ["DEL", "LOSS", "1", "0"]):
        return "LOSS"
    if any(x in t for x in ["DUP", "GAIN", "3", "4"]):
        return "GAIN"
    return "OTHER"

def main():
    args = parse_args()
    
    # 1. Chargement de la verite terrain
    truth_df = pd.read_csv(args.truth, sep=None, engine='python')
    truth_df.columns = [c.lower() for c in truth_df.columns]
    truth_df['chr'] = truth_df['chr'].apply(normalize_chr)
    truth_df['norm_type'] = truth_df['type'].apply(normalize_type)
    
    # 2. Chargement des appels du caller
    calls_df = pd.read_csv(args.calls, sep=None, engine='python')
    calls_df.columns = [c.lower() for c in calls_df.columns]
    calls_df['chr'] = calls_df['chr'].apply(normalize_chr)
    calls_df['norm_type'] = calls_df['type'].apply(normalize_type)
    
    tp, fn, fp = 0, 0, 0
    detected_genes = []
    missed_genes = []
    
    # 3. Evaluation pour chaque gene attendu
    for _, row in truth_df.iterrows():
        g_chr = row['chr']
        g_start = int(row['start'])
        g_end = int(row['end'])
        g_type = row['norm_type']
        g_name = row.get('gene', f"{g_chr}:{g_start}-{g_end}")
        g_len = max(1, g_end - g_start)
        
        # Filtrer les segments candidats sur le meme chromosome et meme sens (GAIN/LOSS)
        candidates = calls_df[(calls_df['chr'] == g_chr) & (calls_df['norm_type'] == g_type)]
        
        gene_detected = False
        for _, c_row in candidates.iterrows():
            c_start = int(c_row['start'])
            c_end = int(c_row['end'])
            
            # Calcul du chevauchement
            overlap_start = max(g_start, c_start)
            overlap_end = min(g_end, c_end)
            overlap_len = max(0, overlap_end - overlap_start)
            
            # Condition de recouvrement >= 50% du gene
            if (overlap_len / g_len) >= args.overlap:
                gene_detected = True
                break
                
        if gene_detected:
            tp += 1
            detected_genes.append(g_name)
        else:
            fn += 1
            missed_genes.append(g_name)
            
    # Calcul de la Sensibilite (Recall)
    sensitivity = (tp / (tp + fn)) * 100 if (tp + fn) > 0 else 0.0
    
    # 4. Ecriture du rapport
    report = f"""================================================================================
BiSHOP V2 — Bilan d'Evaluation CNV (Seuil recouvrement gene: {args.overlap*100:.0f}%)
================================================================================
Genes attendus (Verite)  : {len(truth_df)}
Genes detectes (TP)      : {tp}
Genes manques  (FN)      : {fn}
Sensibilite / Recall     : {sensitivity:.2f} %
================================================================================
Genes detectes : {', '.join(detected_genes) if detected_genes else 'Aucun'}
Genes manques  : {', '.join(missed_genes) if missed_genes else 'Aucun'}
================================================================================
"""
    print(report)
    with open(args.out, "w") as f:
        f.write(report)

if __name__ == "__main__":
    main()
