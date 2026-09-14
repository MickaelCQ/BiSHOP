#!/usr/bin/env python3
"""
Consensus Intelligent SNV/Indel avec Traçabilité des Callers (BiSHOP V2)
Injecte CALLERS=... et NUM_CALLERS=... dans le VCF
"""
import sys
import gzip
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dv", required=True, help="VCF DeepVariant normalisé")
    parser.add_argument("--gatk", required=True, help="VCF GATK normalisé")
    parser.add_argument("--fb", required=True, help="VCF FreeBayes normalisé")
    parser.add_argument("--sample", required=True, help="Nom de l'échantillon")
    parser.add_argument("--out", required=True, help="VCF consensus de sortie")
    return parser.parse_args()

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def main():
    args = parse_args()
    
    variants = defaultdict(lambda: {"callers": set(), "line_info": None, "format": None, "sample": None, "qual": 0.0})
    header_lines = []
    
    callers = [("DeepVariant", args.dv), ("GATK", args.gatk), ("FreeBayes", args.fb)]
    
    for caller_name, vcf_file in callers:
        with open_vcf(vcf_file) as f:
            for line in f:
                if line.startswith("##"):
                    if caller_name == "DeepVariant" and not line.startswith("##INFO=<ID=CALLERS") and not line.startswith("##INFO=<ID=NUM_CALLERS"):
                        header_lines.append(line.strip())
                    continue
                if line.startswith("#CHROM"):
                    continue
                    
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                chrom, pos, var_id, ref, alt, qual, flt, info, fmt, smp = parts[:10]
                key = (chrom, pos, ref, alt)
                
                try:
                    q_val = float(qual) if qual != "." else 0.0
                except ValueError:
                    q_val = 0.0
                    
                variants[key]["callers"].add(caller_name)
                # On privilégie les métriques DeepVariant ou le caller avec le meilleur QUAL
                if variants[key]["line_info"] is None or caller_name == "DeepVariant" or q_val > variants[key]["qual"]:
                    variants[key]["line_info"] = (chrom, pos, var_id, ref, alt, qual, flt, info)
                    variants[key]["format"] = fmt
                    variants[key]["sample"] = smp
                    variants[key]["qual"] = q_val

    # Écriture du VCF Consensus enrichi
    with open(args.out, "w") as out:
        for h in header_lines:
            out.write(h + "\n")
        out.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Callers ayant valide le variant (DeepVariant, GATK, FreeBayes)">\n')
        out.write('##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description="Nombre de callers en accord (1, 2 ou 3)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
        
        # Tri génomique
        sorted_keys = sorted(variants.keys(), key=lambda x: (x[0], int(x[1])))
        for key in sorted_keys:
            data = variants[key]
            chrom, pos, var_id, ref, alt, qual, flt, info = data["line_info"]
            caller_list = ",".join(sorted(data["callers"]))
            num_callers = len(data["callers"])
            
            new_info = f"{info};CALLERS={caller_list};NUM_CALLERS={num_callers}"
            out.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{new_info}\t{data['format']}\t{data['sample']}\n")

if __name__ == "__main__":
    main()
