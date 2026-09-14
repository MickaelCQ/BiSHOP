#!/usr/bin/env python3
"""
Générateur de la Matrice Globale de Concordance Clinique (Multi-Patients)
Panel Custom vs Exome WES — CHU de Nîmes (Norme ISO 15189)
"""
import os
import glob
import pandas as pd

def build_master_matrix(annot_bed, panel_type, out_csv):
    print("=" * 85)
    print(f"GÉNÉRATION DU GRAND TABLEAU DE CONCORDANCE — PANEL : {panel_type}")
    print("=" * 85)
    
    # 1. Charger le BED annoté officiel
    bed_df = pd.read_csv(annot_bed, sep="\t", header=None)
    bed_df = bed_df.iloc[:, [0, 1, 2, 3, 4, 5]]
    bed_df.columns = ["chr", "start", "end", "region_id", "gene", "type"]
    n_regions = len(bed_df)
    
    # 2. Mapping des patients (Legacy vs WES)
    if panel_type == "SLA":
        pairs = [
            ("NA24385-2508121471_SLA", "ADN_CORIELL_NA24385-A", "Coriell_A"),
            ("NA24385-2508121471_SLA", "ADN_CORIELL_NA24385-B", "Coriell_B"),
            ("2107192220_SLA", "2107192220-A", "P_2107192220_A"),
            ("2107192220_SLA", "2107192220-B", "P_2107192220_B"),
            ("2207131857_SLA", "2207131857", "P_2207131857"),
            ("2408070565_SLA", "2408070565", "P_2408070565"),
            ("2409131379_SLA", "2409131379", "P_2409131379"),
            ("2603051686_SLA", "2603051686-60X", "P_2603051686_60X"),
            ("2603051686_SLA", "2603051686-120X", "P_2603051686_120X"),
            ("2603051686_SLA", "2603051686-240X", "P_2603051686_240X")
        ]
    else: # TSA
        pairs = [
            ("2006081825_TSA-HF", "2006081825", "P_2006081825"),
            ("2407011622_TSA-HF", "2407011622", "P_2407011622")
        ]
        
    all_wes_cols = []
    all_leg_cols = []
    
    for leg_id, wes_id, label in pairs:
        leg_file = os.path.expanduser(f"~/Exome/BiSHOP_V2/legacy_panels/qc/{leg_id}_legacy.regions.bed.gz")
        wes_file = os.path.expanduser(f"~/Exome/BiSHOP_V2/legacy_panels/qc/wes/{wes_id}_wes.regions.bed.gz")
        
        if os.path.exists(leg_file) and os.path.exists(wes_file):
            leg_s = pd.read_csv(leg_file, sep="\t", compression="gzip", header=None).iloc[:, -1].values[:n_regions]
            wes_s = pd.read_csv(wes_file, sep="\t", compression="gzip", header=None).iloc[:, -1].values[:n_regions]
            
            c_leg = f"cov_panel_{label}"
            c_wes = f"cov_wes_{label}"
            c_stat = f"statut_{label}"
            
            bed_df[c_leg] = pd.Series(leg_s).round(1)
            bed_df[c_wes] = pd.Series(wes_s).round(1)
            
            # Statut individuel par patient
            def eval_sample(row, cl=c_leg, cw=c_wes):
                if row[cw] >= 100:
                    return "OPTIMAL (>=100X)"
                elif row[cw] >= 30:
                    return "CONFORME (>=30X)"
                elif row[cl] >= 30 and row[cw] < 30:
                    return "DROPOUT_CRITIQUE (<30X)"
                else:
                    return "SOUS_COUVERTURE (<30X)"
                    
            bed_df[c_stat] = bed_df.apply(eval_sample, axis=1)
            all_leg_cols.append(c_leg)
            all_wes_cols.append(c_wes)
            print(f"  [OK] Confrontation intégrée pour : {label}")
        else:
            print(f"  [SKIP] Fichier manquant pour {label} (leg={os.path.exists(leg_file)}, wes={os.path.exists(wes_file)})")
            
    # 3. Calcul des métriques de synthèse de cohorte
    bed_df["Moyenne_Panel_Cohorte"] = bed_df[all_leg_cols].mean(axis=1).round(1)
    bed_df["Moyenne_WES_Cohorte"]   = bed_df[all_wes_cols].mean(axis=1).round(1)
    bed_df["Ratio_Global_WES_Panel"] = (bed_df["Moyenne_WES_Cohorte"] / bed_df["Moyenne_Panel_Cohorte"].replace(0, 0.1)).round(2)
    
    # Statut global d'accréditation pour le gène/exon
    def eval_global(row):
        # Vérifie si le moindre échantillon a un dropout
        stat_cols = [c for c in row.index if c.startswith("statut_")]
        statuses = [row[c] for c in stat_cols]
        if any("DROPOUT_CRITIQUE" in s for s in statuses):
            return "ALERTE_DROPOUT"
        elif row["Moyenne_WES_Cohorte"] >= 100:
            return "VALIDÉ_OPTIMAL (>=100X)"
        elif row["Moyenne_WES_Cohorte"] >= 30:
            return "VALIDÉ_CONFORME (>=30X)"
        else:
            return "ATTENTION_SOUS_COUVERTURE"
            
    bed_df["Synthese_Accreditation"] = bed_df.apply(eval_global, axis=1)
    
    # Réorganisation des colonnes : Infos région -> Synthèse globale -> Détail par patient
    meta_cols = ["gene", "region_id", "chr", "start", "end", "type", 
                 "Synthese_Accreditation", "Moyenne_WES_Cohorte", "Moyenne_Panel_Cohorte", "Ratio_Global_WES_Panel"]
    patient_cols = [c for c in bed_df.columns if c not in meta_cols]
    final_df = bed_df[meta_cols + patient_cols].sort_values(by=["gene", "region_id"])
    
    final_df.to_csv(out_csv, sep=";", index=False)
    
    # Résumé clinique affiché à l'écran
    n_total = len(final_df)
    n_optimal = len(final_df[final_df["Synthese_Accreditation"].str.startswith("VALIDÉ_OPTIMAL")])
    n_conforme = len(final_df[final_df["Synthese_Accreditation"].str.startswith("VALIDÉ_CONFORME")])
    n_dropouts = len(final_df[final_df["Synthese_Accreditation"] == "ALERTE_DROPOUT"])
    
    print("-" * 85)
    print(f"RÉSUMÉ ACCRÉDITATION ISO 15189 — PANEL {panel_type} :")
    print(f"Régions totales diagnostiques        : {n_total}")
    print(f"Exons validés optimaux (>= 100X WES) : {n_optimal} ({n_optimal/n_total*100:.1f} %)")
    print(f"Exons validés conformes (>= 30X WES) : {n_optimal + n_conforme} ({(n_optimal+n_conforme)/n_total*100:.1f} %)")
    print(f"Exons avec alertes Dropouts (< 30X)  : {n_dropouts} ({n_dropouts/n_total*100:.1f} %)")
    print(f"Profondeur moyenne globale WES       : {final_df['Moyenne_WES_Cohorte'].mean():.1f} X")
    print(f"Profondeur moyenne globale Panel     : {final_df['Moyenne_Panel_Cohorte'].mean():.1f} X")
    print(f"Fichier CSV Maître généré            : {out_csv}")
    print("=" * 85 + "\n")

# Lancement pour SLA
build_master_matrix(
    os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels/bed/SLA_Panel_hg38_annotated.bed"),
    "SLA",
    os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels/concordance/SLA_Master_Clinical_Concordance.csv")
)

# Lancement pour TSA
build_master_matrix(
    os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels/bed/TSA_Panel_hg38_annotated.bed"),
    "TSA",
    os.path.expanduser("~/Exome/BiSHOP_V2/legacy_panels/concordance/TSA_Master_Clinical_Concordance.csv")
)
