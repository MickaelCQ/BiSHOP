#!/usr/bin/env python3
"""
Générateur des Deux Fichiers Maîtres de Synthèse Globale (ToBeConcluded)
BiSHOP V2 - CHU de Nîmes (Norme ISO 15189)
"""
import os
import sys
import glob
import gzip
import pandas as pd

def consolidate_snv(snv_csv_files, out_snv_master):
    print("=" * 85)
    print(f">>> 1. CONSOLIDATION GLOBALE SNV / INDEL ({len(snv_csv_files)} fichiers patients)")
    dfs = []
    for f in snv_csv_files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            try:
                d = pd.read_csv(f, sep=";", encoding="utf-8-sig")
                dfs.append(d)
            except Exception as e:
                print(f"  [WARN] Erreur lecture {f}: {e}")
                
    if dfs:
        master_snv = pd.concat(dfs, ignore_index=True)
        # Tri par patient, puis chromosome et position
        master_snv["Pos_num"] = pd.to_numeric(master_snv["Pos"], errors="coerce").fillna(0).astype(int)
        master_snv = master_snv.sort_values(by=["Sample", "Chr", "Pos_num"]).drop(columns=["Pos_num"])
        master_snv.to_csv(out_snv_master, sep=";", index=False, encoding="utf-8-sig")
        print(f"  [OK] Master SNV/INDEL généré : {len(master_snv)} variants répertoriés dans {out_snv_master}")
    else:
        print("  [WARN] Aucun fichier SNV trouvé.")

def consolidate_cnv(clincnv_files, decon_csv, cnvkit_files, out_cnv_master):
    print("=" * 85)
    print(">>> 2. CONSOLIDATION GLOBALE DES CNVs (CNVkit + ClinCNV + DECoN)")
    cnv_records = []

    # A. Récupérer ClinCNV (TSV de chaque patient)
    for f in clincnv_files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            sname = os.path.basename(f).replace("_cnvs.tsv", "")
            try:
                df = pd.read_csv(f, sep="\t", comment="#")
                for _, r in df.iterrows():
                    c = str(r.get("chr", r.get("#chr", ""))).replace("chr", "")
                    s, e = int(r["start"]), int(r["end"])
                    cn = int(r.get("CN_change", 2))
                    svtype = "DEL" if cn < 2 else "DUP" if cn > 2 else "NORMAL"
                    genes = str(r.get("genes", "."))
                    qv = r.get("qvalue", ".")
                    if svtype != "NORMAL":
                        cnv_records.append({
                            "Sample": sname, "Caller": "ClinCNV", "Chr": c,
                            "Start": s, "End": e, "Size_kb": round((e - s) / 1000, 2),
                            "Type": svtype, "Copy_Number": cn, "Genes": genes,
                            "Score_Confiance": f"qval={qv}"
                        })
            except Exception: pass

    # B. Récupérer DECoN
    if os.path.exists(decon_csv) and os.path.getsize(decon_csv) > 0:
        try:
            df = pd.read_csv(decon_csv, sep=",")
            for _, r in df.iterrows():
                raw_s = str(r.get("sample.names.i.", r.get("Sample", "")))
                sname = raw_s.replace(".markdup.bam", "").replace(".sorted.bam", "").replace(".bam", "")
                c = str(r.get("chromosome", r.get("Chromosome", ""))).replace("chr", "")
                s, e = int(r.get("start", r.get("Start", 0))), int(r.get("end", r.get("End", 0)))
                t = str(r.get("type", r.get("Type", ""))).upper()
                svtype = "DEL" if "DEL" in t else "DUP" if "DUP" in t else "OTHER"
                genes = str(r.get("name", r.get("Exons", ".")))
                bf = r.get("BF", ".")
                if svtype != "OTHER":
                    cnv_records.append({
                        "Sample": sname, "Caller": "DECoN", "Chr": c,
                        "Start": s, "End": e, "Size_kb": round((e - s) / 1000, 2),
                        "Type": svtype, "Copy_Number": "1 (hét)" if svtype == "DEL" else "3 (dup)",
                        "Genes": genes, "Score_Confiance": f"BF={bf}"
                    })
        except Exception: pass

    # C. Récupérer CNVkit (VCFs)
    for f in cnvkit_files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            sname = os.path.basename(f).replace(".cnvkit.vcf.gz", "").replace(".cnvkit.vcf", "")
            fh = gzip.open(f, "rt") if f.endswith(".gz") else open(f, "r")
            for line in fh:
                if line.startswith("#"): continue
                p = line.strip().split("\t")
                c = p[0].replace("chr", "")
                s = int(p[1])
                e = s
                svtype = "DEL" if "<DEL>" in p[4] else "DUP" if "<DUP>" in p[4] else "OTHER"
                cn = "."
                for it in p[7].split(";"):
                    if it.startswith("END="): e = int(it.split("=")[1])
                    if it.startswith("CN="): cn = it.split("=")[1]
                if svtype != "OTHER":
                    cnv_records.append({
                        "Sample": sname, "Caller": "CNVkit", "Chr": c,
                        "Start": s, "End": e, "Size_kb": round((e - s) / 1000, 2),
                        "Type": svtype, "Copy_Number": cn, "Genes": ".",
                        "Score_Confiance": f"QUAL={p[5]}"
                    })
            fh.close()

    if cnv_records:
        df_cnv = pd.DataFrame(cnv_records)
        df_cnv = df_cnv.sort_values(by=["Sample", "Chr", "Start"])
        df_cnv["Lien_Decipher"] = "https://www.deciphergenomics.org/browser#q/" + df_cnv["Chr"] + ":" + df_cnv["Start"].astype(str) + "-" + df_cnv["End"].astype(str)
        df_cnv.to_csv(out_cnv_master, sep=";", index=False, encoding="utf-8-sig")
        print(f"  [OK] Master CNV généré : {len(df_cnv)} événements répertoriés dans {out_cnv_master}")
    else:
        print("  [WARN] Aucun événement CNV trouvé.")

def main():
    if len(sys.argv) < 3:
        print("Usage: python consolidate_run_tables.py <results_dir> <out_dir>")
        sys.exit(1)
        
    res_dir = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    
    out_snv = os.path.join(out_dir, "Cohort_Master_SNV_INDEL.csv")
    out_cnv = os.path.join(out_dir, "Cohort_Master_CNV.csv")
    
    snv_files = sorted(glob.glob(os.path.join(res_dir, "annotated_reports", "*.clinical_variants.csv")))
    clincnv_files = sorted(glob.glob(os.path.join(res_dir, "vcfs", "clincnv", "*_cnvs.tsv")))
    decon_csv = os.path.join(res_dir, "vcfs", "decon", "cohort.decon.csv")
    cnvkit_files = sorted(glob.glob(os.path.join(res_dir, "vcfs", "cnvkit", "*.cnvkit.vcf.gz")))
    
    consolidate_snv(snv_files, out_snv)
    consolidate_cnv(clincnv_files, decon_csv, cnvkit_files, out_cnv)
    print("=" * 85)
    print(f"DOSSIER D'ANALYSE GLOBALE GÉNÉRÉ DANS : {out_dir}")
    print("=" * 85 + "\n")

if __name__ == "__main__":
    main()
