#!/usr/bin/env python3
import gzip
from pyliftover import LiftOver

lo = LiftOver('legacy_panels/hg19ToHg38.over.chain.gz')

def lift_bed(input_bed, output_bed):
    converted, unmapped = 0, 0
    with open(input_bed, 'r') as f_in, open(output_bed, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            rest = parts[3:] if len(parts) > 3 else []
            
            c_name = chrom if chrom.startswith('chr') else 'chr' + chrom
            new_start = lo.convert_coordinate(c_name, start)
            new_end = lo.convert_coordinate(c_name, end)
            
            if new_start and new_end and new_start[0][0] == new_end[0][0]:
                out_chr = new_start[0][0].replace('chr', '')
                out_start = new_start[0][1]
                out_end = new_end[0][1]
                f_out.write(f"{out_chr}\t{out_start}\t{out_end}" + ("\t" + "\t".join(rest) if rest else "") + "\n")
                converted += 1
            else:
                unmapped += 1
    print(f"[{input_bed}] -> {converted} régions converties en hg38, {unmapped} non mappées.")

lift_bed('legacy_panels/bed/SLA_Panel.bed', 'legacy_panels/bed/SLA_Panel_hg38.bed')
lift_bed('legacy_panels/bed/TSA_Panel.bed', 'legacy_panels/bed/TSA_Panel_hg38.bed')
