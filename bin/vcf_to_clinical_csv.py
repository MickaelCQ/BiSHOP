#!/usr/bin/env python3
"""
Convertisseur Officiel VCF Annoté VEP -> CSV Médical Clinique
BiSHOP V2 - CHU de Nîmes (Norme ISO 15189)
"""
import sys
import glob
import gzip
import urllib.parse
import pandas as pd

def open_file(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f, "r")

def main():
    if len(sys.argv) < 3:
        print("Usage: python vcf_to_clinical_csv.py <input.vep.vcf.gz> <output.csv> [sample_id] [activity]")
        sys.exit(1)

    raw_path = sys.argv[1]
    matched = glob.glob(raw_path)
    if not matched:
        print(f"Erreur : Aucun fichier trouvé correspondant à '{raw_path}'")
        sys.exit(1)
    vcf_path = matched[0]

    out_csv = sys.argv[2]
    sample_id = sys.argv[3] if len(sys.argv) > 3 else "Sample"
    activity = sys.argv[4] if len(sys.argv) > 4 else "WES"

    csq_fields = []
    with open_file(vcf_path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                desc = line.split("Format: ")[1].split('">')[0]
                csq_fields = desc.split("|")
                break
            if not line.startswith("#"): break

    csq_map = {name: idx for idx, name in enumerate(csq_fields)}
    rows = []

    with open_file(vcf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 10: continue

            chrom, pos, var_id, ref, alt, qual, flt, info_str, fmt_str, smp_str = parts[:10]

            # Extraction INFO
            info_dict = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            callers = info_dict.get("CALLERS", ".")
            num_callers = info_dict.get("NUM_CALLERS", ".")

            # Extraction FORMAT & SAMPLE
            fmt_keys = fmt_str.split(":")
            smp_vals = smp_str.split(":")
            smp_dict = dict(zip(fmt_keys, smp_vals))

            gt = smp_dict.get("GT", ".")
            dp = smp_dict.get("DP", ".")
            ad = smp_dict.get("AD", ".")
            vaf = smp_dict.get("VAF", ".")

            if (vaf == "." or vaf == "") and ad != "." and "," in ad:
                try:
                    ref_c, alt_c = [float(x) for x in ad.split(",")[:2]]
                    tot = ref_c + alt_c
                    if tot > 0: vaf = f"{alt_c / tot:.3f}"
                except ValueError: pass

            # Extraction VEP CSQ (Transcrit Canonique prioritaire)
            csq_str = info_dict.get("CSQ", "")
            transcripts = csq_str.split(",") if csq_str else []
            selected_f = None
            
            for tr in transcripts:
                fields = tr.split("|")
                if len(fields) < len(csq_fields):
                    fields += [""] * (len(csq_fields) - len(fields))
                is_canonical = fields[csq_map.get("CANONICAL", -1)] if "CANONICAL" in csq_map else ""
                if is_canonical == "YES":
                    selected_f = fields
                    break

            if not selected_f and transcripts:
                selected_f = transcripts[0].split("|")
                if len(selected_f) < len(csq_fields):
                    selected_f += [""] * (len(csq_fields) - len(selected_f))

            def get_csq(tag):
                if not selected_f or tag not in csq_map: return "."
                idx = csq_map[tag]
                val = selected_f[idx] if idx < len(selected_f) else ""
                return urllib.parse.unquote(val) if val else "."

            gene = get_csq("SYMBOL")
            biotype = get_csq("BIOTYPE")
            tr_id = get_csq("Feature")
            hgvsc = get_csq("HGVSc")
            hgvsp = get_csq("HGVSp")
            csq = get_csq("Consequence")
            impact = get_csq("IMPACT")
            exon = get_csq("EXON") if get_csq("EXON") != "." else get_csq("INTRON")
            
            # Extraction des bases cliniques officielles
            clinvar = get_csq("ClinVar_CLNSIG") if get_csq("ClinVar_CLNSIG") != "." else info_dict.get("ClinVar_CLNSIG", get_csq("CLIN_SIG"))
            clin_rev = get_csq("ClinVar_CLNREVSTAT") if get_csq("ClinVar_CLNREVSTAT") != "." else info_dict.get("ClinVar_CLNREVSTAT", ".")
            am_score = get_csq("AlphaMissense_am_pathogenicity") if get_csq("AlphaMissense_am_pathogenicity") != "." else info_dict.get("AlphaMissense_am_pathogenicity", ".")
            am_class = get_csq("AlphaMissense_am_class") if get_csq("AlphaMissense_am_class") != "." else info_dict.get("AlphaMissense_am_class", ".")
            revel = get_csq("REVEL")
            gnomad = get_csq("gnomADe_AF") if get_csq("gnomADe_AF") != "." else get_csq("AF")

            # Lien direct MobiDetails
            mobidetails_url = f"https://mobidetails.iurc.montp.inserm.fr/MD/variant/?variant=chr{chrom}-{pos}-{ref}-{alt}&genome_version=hg38"

            rows.append({
                "Sample": sample_id,
                "Activite": activity,
                "Chr": chrom,
                "Pos": pos,
                "Ref": ref,
                "Alt": alt,
                "Qual": qual,
                "Filter": flt,
                "Callers": callers,
                "Nb_Callers": num_callers,
                "Genotype": gt,
                "Profondeur_DP": dp,
                "Reads_AD": ad,
                "VAF": vaf,
                "Gene": gene,
                "Biotype": biotype,
                "Transcript_ID": tr_id,
                "HGVSc": hgvsc,
                "HGVSp": hgvsp,
                "Consequence": csq,
                "Impact": impact,
                "Exon_Intron": exon,
                "ClinVar_Significance": clinvar,
                "ClinVar_ReviewStatus": clin_rev,
                "AlphaMissense_Score": am_score,
                "AlphaMissense_Class": am_class,
                "REVEL_Score": revel,
                "gnomAD_AF": gnomad,
                "Lien_MobiDetails": mobidetails_url
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, sep=";", index=False, encoding="utf-8-sig")
    print(f"[{sample_id}] Tableau clinique généré avec succès : {len(df)} variants dans {out_csv}")

if __name__ == "__main__":
    main()
