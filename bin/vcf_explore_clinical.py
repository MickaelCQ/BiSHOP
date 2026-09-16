#!/usr/bin/env python3
"""
Explorateur Clinique VCF & VEP Dynamique (Compatible VisiData & Excel)
Génère des colonnes explicites avec les options entre parenthèses
BiSHOP V2 - CHU de Nîmes
"""
import sys
import gzip
import urllib.parse

def open_file(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f, "r")

def main():
    if len(sys.argv) < 2:
        print("Usage: python vcf_explore_clinical.py <fichier.vep.vcf.gz> [target_gene]")
        sys.exit(1)

    vcf_path = sys.argv[1]
    target_gene = sys.argv[2] if len(sys.argv) > 2 else ""

    csq_fields = []
    
    # 1. Analyse dynamique de l'en-tête pour cartographier les barres | de VEP
    with open_file(vcf_path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                desc = line.split("Format: ")[1].split('">')[0]
                csq_fields = desc.split("|")
                break
            if not line.startswith("#"):
                break

    csq_map = {name: idx for idx, name in enumerate(csq_fields)}

    # 2. Définition des en-têtes demandés avec les options entre parenthèses
    headers = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "CALLERS (--consensus)",
        "NUM_CALLERS (--consensus)",
        "GT (Genotype)",
        "DP (Profondeur)",
        "AD (Reads_Ref,Alt)",
        "VAF (Frequence_Allel)",
        "GENE (--symbol)",
        "BIOTYPE (--biotype)",
        "TRANSCRIPT (--canonical)",
        "HGVSc (--hgvs)",
        "HGVSp (--hgvs)",
        "CONSEQUENCE (--everything)",
        "IMPACT (--everything)",
        "EXON/INTRON (--numbers)",
        "CLIN_SIG (--clin_sig_allele)",
        "REVEL (--plugin REVEL)",
        "GNOMAD_AF (--af_gnomad)"
    ]
    print("\t".join(headers))

    # 3. Traitement des variants
    with open_file(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

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

            # Calcul dynamique du VAF si manquant (notamment pour FreeBayes)
            if (vaf == "." or vaf == "") and ad != "." and "," in ad:
                try:
                    ref_c, alt_c = [float(x) for x in ad.split(",")[:2]]
                    tot = ref_c + alt_c
                    if tot > 0:
                        vaf = f"{alt_c / tot:.4f}"
                except ValueError:
                    pass

            # Extraction dynamique VEP dans le champ CSQ
            csq_str = info_dict.get("CSQ", "")
            transcripts = csq_str.split(",") if csq_str else []

            # Recherche du transcrit canonique ou correspondant au gène cible
            selected_f = None
            for tr in transcripts:
                fields = tr.split("|")
                # Compléter avec des chaînes vides si besoin
                if len(fields) < len(csq_fields):
                    fields += [""] * (len(csq_fields) - len(fields))

                gene_name = fields[csq_map.get("SYMBOL", 3)]
                is_canonical = fields[csq_map.get("CANONICAL", -1)] if "CANONICAL" in csq_map else ""

                if target_gene and gene_name.upper() == target_gene.upper():
                    selected_f = fields
                    break
                elif is_canonical == "YES":
                    selected_f = fields
                    break

            if not selected_f and transcripts:
                selected_f = transcripts[0].split("|")
                if len(selected_f) < len(csq_fields):
                    selected_f += [""] * (len(csq_fields) - len(selected_f))

            # Récupération des valeurs annotées
            def get_csq(tag):
                if not selected_f or tag not in csq_map:
                    return "."
                idx = csq_map[tag]
                val = selected_f[idx] if idx < len(selected_f) else ""
                val = urllib.parse.unquote(val) # Décode %3D en =
                return val if val else "."

            gene = get_csq("SYMBOL")
            biotype = get_csq("BIOTYPE")
            tr_id = get_csq("Feature")
            hgvsc = get_csq("HGVSc")
            hgvsp = get_csq("HGVSp")
            csq = get_csq("Consequence")
            impact = get_csq("IMPACT")
            exon = get_csq("EXON") if get_csq("EXON") != "." else get_csq("INTRON")
            clin_sig = get_csq("CLIN_SIG")
            revel = get_csq("REVEL")
            gnomad = get_csq("gnomADe_AF") if get_csq("gnomADe_AF") != "." else get_csq("AF")

            if target_gene and gene.upper() != target_gene.upper():
                continue

            row = [
                chrom, pos, ref, alt, qual, flt,
                callers, num_callers,
                gt, dp, ad, vaf,
                gene, biotype, tr_id,
                hgvsc, hgvsp, csq, impact,
                exon, clin_sig, revel, gnomad
            ]
            print("\t".join(row))

if __name__ == "__main__":
    main()
