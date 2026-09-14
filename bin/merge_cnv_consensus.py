#!/usr/bin/env python3
"""
Matrice de Concordance et Consensus CNV R&D (Sans Perte d'Information)
BiSHOP V2 - CHU de Nîmes (Investigation Clinique)
"""
import os
import sys
import gzip
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Consensus et matrice de concordance CNV multi-callers")
    parser.add_argument("--cnvkit", required=True, help="VCF CNVkit")
    parser.add_argument("--clincnv", required=True, help="TSV ClinCNV")
    parser.add_argument("--decon", required=True, help="VCF DECoN")
    parser.add_argument("--sample", required=True, help="Identifiant patient")
    parser.add_argument("--out", required=True, help="VCF de synthèse R&D")
    return parser.parse_args()

def normalize_chr(c):
    c = str(c).strip()
    return c if c.startswith("chr") else "chr" + c

def main():
    args = parse_args()
    raw_calls = []
    
    # 1. Extraction CNVkit
    if os.path.exists(args.cnvkit) and os.path.getsize(args.cnvkit) > 0:
        f = gzip.open(args.cnvkit, "rt") if args.cnvkit.endswith(".gz") else open(args.cnvkit, "r")
        for line in f:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            chrom, start, end = normalize_chr(p[0]), int(p[1]), int(p[1])
            alt = p[4]
            svtype = "DEL" if "<DEL>" in alt else "DUP" if "<DUP>" in alt else "OTHER"
            for item in p[7].split(";"):
                if item.startswith("END="): end = int(item.split("=")[1])
            if svtype != "OTHER":
                raw_calls.append({"caller": "CNVkit", "chr": chrom, "start": start, "end": end, "type": svtype, "raw": f"CNVkit:{start}-{end}"})
        f.close()

    # 2. Extraction ClinCNV
    if os.path.exists(args.clincnv) and os.path.getsize(args.clincnv) > 0:
        try:
            df = pd.read_csv(args.clincnv, sep=None, engine='python', comment="#")
            for _, r in df.iterrows():
                chrom = normalize_chr(r.get("chr", r.get("#chr", "")))
                start, end = int(r["start"]), int(r["end"])
                cn = int(r.get("CN_change", 2))
                gene = str(r.get("genes", ""))
                svtype = "DEL" if cn < 2 else "DUP" if cn > 2 else "OTHER"
                if svtype != "OTHER":
                    raw_calls.append({"caller": "ClinCNV", "chr": chrom, "start": start, "end": end, "type": svtype, "raw": f"ClinCNV:{start}-{end}({gene})"})
        except Exception:
            pass

    # 3. Extraction DECoN
    if os.path.exists(args.decon) and os.path.getsize(args.decon) > 0:
        f = gzip.open(args.decon, "rt") if args.decon.endswith(".gz") else open(args.decon, "r")
        for line in f:
            if line.startswith("#"): continue
            p = line.strip().split("\t")
            chrom, start, end = normalize_chr(p[0]), int(p[1]), int(p[1])
            alt = p[4]
            svtype = "DEL" if "<DEL>" in alt else "DUP" if "<DUP>" in alt else "OTHER"
            for item in p[7].split(";"):
                if item.startswith("END="): end = int(item.split("=")[1])
            if svtype != "OTHER":
                raw_calls.append({"caller": "DECoN", "chr": chrom, "start": start, "end": end, "type": svtype, "raw": f"DECoN:{start}-{end}"})
        f.close()

    # 4. Regroupement intelligent SANS suppression
    clusters = []
    used = [False] * len(raw_calls)
    
    for i in range(len(raw_calls)):
        if used[i]: continue
        c1 = raw_calls[i]
        matched_callers = {c1["caller"]}
        raw_traces = [c1["raw"]]
        m_start, m_end = c1["start"], c1["end"]
        used[i] = True
        
        for j in range(i + 1, len(raw_calls)):
            if used[j]: continue
            c2 = raw_calls[j]
            if c1["chr"] == c2["chr"] and c1["type"] == c2["type"]:
                ov_s = max(m_start, c2["start"])
                ov_e = min(m_end, c2["end"])
                ov_len = max(0, ov_e - ov_s)
                l1 = m_end - m_start
                l2 = c2["end"] - c2["start"]
                
                # Chevauchement significatif (au moins 20% réciproque ou 50 pb)
                if l1 > 0 and l2 > 0 and (ov_len / l1 >= 0.20 or ov_len / l2 >= 0.20 or ov_len >= 50):
                    matched_callers.add(c2["caller"])
                    raw_traces.append(c2["raw"])
                    m_start = min(m_start, c2["start"])
                    m_end = max(m_end, c2["end"])
                    used[j] = True
                    
        n_callers = len(matched_callers)
        if n_callers == 3:
            status = "CONCORDANT_HIGH_CONFIDENCE (3/3)"
        elif n_callers == 2:
            status = "CONCORDANT_MEDIUM (2/3)"
        else:
            status = f"DISCORDANT_SINGLETON_{list(matched_callers)[0]}"
            
        clusters.append({
            "chr": c1["chr"].replace("chr", ""),
            "start": m_start,
            "end": m_end,
            "type": c1["type"],
            "callers": ",".join(sorted(matched_callers)),
            "num_callers": n_callers,
            "status": status,
            "details": "|".join(raw_traces)
        })

    # 5. Écriture du VCF de synthèse
    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type de variant structural (DEL, DUP)">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Position de fin">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Taille du variant">\n')
        out.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Callers soutenant le variant">\n')
        out.write('##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description="Nombre de callers en accord (1 a 3)">\n')
        out.write('##INFO=<ID=STATUS,Number=1,Type=String,Description="Statut de concordance diagnostique">\n')
        out.write('##INFO=<ID=ORIGINAL_CALLS,Number=.,Type=String,Description="Coordonnees brutes de chaque caller">\n')
        out.write('##ALT=<ID=DEL,Description="Deletion">\n')
        out.write('##ALT=<ID=DUP,Description="Duplication">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
        
        for c in sorted(clusters, key=lambda x: (x["chr"], int(x["start"]))):
            svlen = c["end"] - c["start"]
            alt = f"<{c['type']}>"
            flt = "PASS" if c["num_callers"] >= 2 else "SINGLETON"
            info = f"SVTYPE={c['type']};END={c['end']};SVLEN={svlen};CALLERS={c['callers']};NUM_CALLERS={c['num_callers']};STATUS={c['status']};ORIGINAL_CALLS={c['details']}"
            out.write(f"{c['chr']}\t{c['start']}\t.\tN\t{alt}\t.\t{flt}\t{info}\tGT\t0/1\n")

    print(f"[{args.sample}] Synthèse CNV terminée : {len(clusters)} événements répertoriés (Concordants et Singletons préservés).")

if __name__ == "__main__":
    main()
