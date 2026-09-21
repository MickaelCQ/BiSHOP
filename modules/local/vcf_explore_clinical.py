#!/usr/bin/env python3
import sys
import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--activity', required=True)
    parser.add_argument('--out', required=True)
    return parser.parse_args()

def get_csq_header(vcf_path):
    # Extraction de l'en-tête dynamique du cache VEP
    with (gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')) as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                desc = line.strip().split('Format: ')[1].rstrip('">')
                return desc.split('|')
    return []

def main():
    args = parse_args()
    csq_fields = get_csq_header(args.vcf)
    
    with open(args.out, 'w') as out_f, (gzip.open(args.vcf, 'rt') if args.vcf.endswith('.gz') else open(args.vcf, 'r')) as in_f:
        
        # En-têtes optimisés pour lecture médicale
        headers = [
            "Sample", "Activity", "Chr", "Pos", "Ref", "Alt", "Qual", "Filter",
            "Callers", "Nb_Callers", "Genotype", "Depth", "Reads_AD", "VAF",
            "Gene", "Biotype", "Transcript_ID", "HGVSc", "HGVSp", 
            "Consequence", "Impact", "Exon_Intron", 
            "ClinVar_Sig", "ClinVar_RevStat", "AlphaMissense_Score", "AlphaMissense_Class", 
            "REVEL_Score", "gnomAD_AF", "MobiDetails_URL"
        ]
        out_f.write('\t'.join(headers) + '\n')
        
        for line in in_f:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filt, info, fmt, sample_dat = cols[:10]
            
            # Parsing du champ INFO
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True
                    
            callers = info_dict.get('CALLERS', '.')
            nb_callers = info_dict.get('NUM_CALLERS', '0')
            
            # Parsing du format Genotype (Empêche le formatage Date Excel)
            fmt_keys = fmt.split(':')
            s_vals = sample_dat.split(':')
            fmt_dict = dict(zip(fmt_keys, s_vals))
            
            gt = fmt_dict.get('GT', './.')
            if gt in ['0/1', '1/0']: gt_str = 'HET (0/1)'
            elif gt in ['1/1']:      gt_str = 'HOM (1/1)'
            elif gt in ['0/0']:      gt_str = 'WT (0/0)'
            else:                    gt_str = gt
            
            dp = fmt_dict.get('DP', '.')
            ad = fmt_dict.get('AD', '.')
            
            # Calcul de la VAF (Variant Allele Frequency)
            vaf = "."
            if ad != "." and "," in ad:
                ads = ad.split(',')
                try:
                    ref_ad, alt_ad = int(ads[0]), int(ads[1])
                    if ref_ad + alt_ad > 0:
                        vaf = str(round(alt_ad / (ref_ad + alt_ad), 4))
                except: pass
            
            # Parsing de VEP (CSQ)
            csq_str = info_dict.get('CSQ', '')
            if not csq_str: continue
                
            # On prend le premier transcrit (le plus pertinent / canonique)
            first_csq = csq_str.split(',')[0].split('|')
            csq = dict(zip(csq_fields, first_csq))
            
            gene = csq.get('SYMBOL', '.')
            biotype = csq.get('BIOTYPE', '.')
            transcript = csq.get('Feature', '.')
            hgvsc = csq.get('HGVSc', '.')
            hgvsp = csq.get('HGVSp', '.')
            consequence = csq.get('Consequence', '.')
            impact = csq.get('IMPACT', '.')
            
            # Empêche Excel de transformer l'exon "10/19" en Date "oct-19"
            exon = csq.get('EXON', '')
            intron = csq.get('INTRON', '')
            if exon: exon_intron = f"Exon {exon}"
            elif intron: exon_intron = f"Intron {intron}"
            else: exon_intron = "."
                
            clinvar_sig = csq.get('CLNSIG', '.') or '.'
            clinvar_rev = csq.get('CLNREVSTAT', '.') or '.'
            am_score = csq.get('am_pathogenicity', '.') or '.'
            am_class = csq.get('am_class', '.') or '.'
            revel = csq.get('REVEL', '.') or '.'
            
            # GnomAD (VEP peut l'écrire sous gnomAD_AF ou gnomADe_AF)
            gnomad = csq.get('gnomAD_AF', '.')
            if gnomad == '.': gnomad = csq.get('gnomADe_AF', '.')
            if gnomad == '.': gnomad = csq.get('gnomADg_AF', '.')
            if not gnomad: gnomad = '.'
            
            # Lien MobiDetails intelligent (basé sur la nomenclature HGVSc si possible)
            if hgvsc != '.' and ':' in hgvsc:
                # Utilise la recherche HGVS officielle de MobiDetails
                clean_hgvsc = hgvsc.split(':')[-1] # Extrait c.1843C>T
                url = f"https://mobidetails.iurc.montp.inserm.fr/MD/variant_search/?variant={transcript}:{clean_hgvsc}"
            else:
                # Variante génomique de secours
                url = f"https://mobidetails.iurc.montp.inserm.fr/MD/variant_search/?variant={chrom}:g.{pos}{ref}%3E{alt}"
            
            row = [
                args.sample, args.activity, chrom, pos, ref, alt, qual, filt,
                callers, nb_callers, gt_str, str(dp), str(ad), str(vaf),
                gene, biotype, transcript, hgvsc, hgvsp, 
                consequence, impact, exon_intron,
                clinvar_sig, clinvar_rev, am_score, am_class, str(revel),
                str(gnomad), url
            ]
            
            out_f.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    main()
