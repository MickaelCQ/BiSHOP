#!/bin/bash
set -euo pipefail

CHAIN="legacy_panels/hg19ToHg38.over.chain.gz"
RAW_ONCO="bin/process_onco_panels_hg19.bed"
DIR_LEGACY="legacy_panels/bam/ONCO"
DIR_WES="results_nextseq_All_from_scratch/aligned_bams"
OUT="legacy_panels/benchmark_onco_technique"
mkdir -p "$OUT"

echo "=== 1. Homogénéisation du BED hg19 ==="
python3 - << 'EOF'
import re

with open("bin/process_onco_panels_hg19.bed", "r") as fin, open("legacy_panels/benchmark_onco_technique/tmp_cleaned.bed", "w") as fout:
    for line in fin:
        line = line.strip()
        if not line or line.startswith(("track", "#", "browser")):
            continue
        p = line.split("\t")
        if len(p) < 3: continue
        c, s, e = p[0], int(p[1]), int(p[2])
        if s == e: e = s + 1
        elif s > e: s, e = e, s
        
        full = "\t".join(p)
        m = re.search(r'GENE_ID=([^;,\t]+)', full)
        if m: gene = m.group(1)
        elif len(p) >= 6 and "_" in p[5]: gene = p[5].split("_")[0]
        elif len(p) >= 4 and "_" in p[3]: gene = p[3].split("_")[0]
        else: gene = p[3] if len(p) >= 4 else "ONCO"
        fout.write(f"{c}\t{s}\t{e}\t{gene}\n")
EOF

# Fusion et création d'un identifiant stable
sort -k1,1 -k2,2n "$OUT/tmp_cleaned.bed" | bedtools merge -i - -c 4 -o distinct | \
awk 'BEGIN{OFS="\t"} {
    print $1, $2, $3, "ONCO_"NR"_" $4, $4
}' > "$OUT/ONCO_Targets_hg19.bed"
rm -f "$OUT/tmp_cleaned.bed"

echo "=== 2. Annotation CDS Canonical RefSeq hg19 ==="
CDS_HG19="legacy_panels/refseq_cds_hg19.bed"
if [ ! -f "$CDS_HG19" ]; then
    echo "Récupération des CDS NCBI RefSeq hg19..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqCurated.txt.gz" | \
    gunzip -c | awk -F'\t' 'BEGIN{OFS="\t"} {
        split($10, starts, ","); split($11, ends, ",");
        for (i=1; i<=$9; i++) {
            s = (starts[i] > $7) ? starts[i] : $7;
            e = (ends[i] < $8) ? ends[i] : $8;
            if (s < e) print $3, s, e, $13;
        }
    }' | sort -k1,1 -k2,2n | bedtools merge -i - > "$CDS_HG19"
fi

# Flag CDS vs HORS_CDS
bedtools coverage -a "$OUT/ONCO_Targets_hg19.bed" -b "$CDS_HG19" | \
awk 'BEGIN{OFS="\t"} {
    cds_status = ($NF >= 0.20) ? "CDS_RefSeq" : "HORS_CDS_SOMATIQUE";
    print $1, $2, $3, $4, $5, cds_status
}' > "$OUT/ONCO_Targets_Annotated_hg19.bed"

echo "=== 3. LiftOver vers hg38 ==="
liftOver "$OUT/ONCO_Targets_Annotated_hg19.bed" "$CHAIN" "$OUT/ONCO_Targets_Annotated_hg38.bed" "$OUT/unmapped.bed"

# Version sans 'chr' impérative pour les BAMs Exome NextSeq 2000
sed 's/^chr//' "$OUT/ONCO_Targets_Annotated_hg38.bed" > "$OUT/ONCO_Targets_Annotated_hg38_noChr.bed"

echo "=== 4. Calcul mosdepth sur les BAMs Legacy (hg19 avec chr) ==="
for b in "$DIR_LEGACY"/*.bam; do
    [ -f "$b" ] || continue
    name=$(basename "$b" .bam)
    echo "Mosdepth Routine : $name ..."
    mosdepth -t 4 -b "$OUT/ONCO_Targets_Annotated_hg19.bed" -n "${OUT}/${name}" "$b"
    zcat "${OUT}/${name}.regions.bed.gz" > "${OUT}/${name}.tsv"
done

echo "=== 5. Calcul mosdepth sur Coriell A & B WES (hg38 sans chr) ==="
for c_bam in "$DIR_WES"/ADN_CORIELL_NA24385-A.markdup.bam "$DIR_WES"/ADN_CORIELL_NA24385-B.markdup.bam; do
    [ -f "$c_bam" ] || continue
    c_name=$(basename "$c_bam" .markdup.bam)
    echo "Mosdepth WES : $c_name ..."
    # On utilise la version noChr pour matcher le BAM Exome
    mosdepth -t 4 -b "$OUT/ONCO_Targets_Annotated_hg38_noChr.bed" -n "${OUT}/${c_name}" "$c_bam"
    zcat "${OUT}/${c_name}.regions.bed.gz" > "${OUT}/${c_name}.tsv"
done

echo "=== 6. Génération du CSV Exploratoire Technique ==="
python3 - << 'EOF'
import os

out = "legacy_panels/benchmark_onco_technique"

# 1. Structure de base depuis le BED hg38 annoté
regions = {}
with open(f"{out}/ONCO_Targets_Annotated_hg38.bed", "r") as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 6: continue
        c, s, e, reg_id, gene, status = p[0], p[1], p[2], p[3], p[4], p[5]
        regions[reg_id] = {
            "chr_hg38": c, "start_hg38": s, "end_hg38": e,
            "gene": gene, "status": status,
            "depths": {}
        }

# 2. Récupération des profondeurs des 4 BAMs panels (dernière colonne p[-1])
legacy_samples = ["2608120547_BRCA", "2608260444_MYELOIDE", "2608280855_OPA-CF", "2608311668_OPA"]
for s in legacy_samples:
    tsv = f"{out}/{s}.tsv"
    if os.path.exists(tsv):
        with open(tsv) as f:
            for l in f:
                p = l.strip().split("\t")
                if len(p) >= 5:
                    reg_id, depth = p[3], float(p[-1])
                    if reg_id in regions:
                        regions[reg_id]["depths"][s] = depth

# 3. Récupération des profondeurs Coriell A et B
coriell_samples = ["ADN_CORIELL_NA24385-A", "ADN_CORIELL_NA24385-B"]
for c in coriell_samples:
    tsv = f"{out}/{c}.tsv"
    if os.path.exists(tsv):
        with open(tsv) as f:
            for l in f:
                p = l.strip().split("\t")
                if len(p) >= 5:
                    reg_id, depth = p[3], float(p[-1])
                    if reg_id in regions:
                        regions[reg_id]["depths"][c] = depth

# 4. Écriture du CSV
csv_file = f"{out}/ONCO_Exploration_Technique_NextSeq.csv"
headers = [
    "gene", "region_id", "chr_hg38", "start_hg38", "end_hg38", "type_region",
    "cov_BRCA_routine", "cov_MYELOIDE_routine", "cov_OPA_CF_routine", "cov_OPA_routine",
    "Max_Baseline_Routine",
    "cov_NextSeq_Coriell_A", "cov_NextSeq_Coriell_B", "Moyenne_NextSeq_Coriell",
    "Ratio_NextSeq_vs_Baseline", "Diagnostic_Technique"
]

with open(csv_file, "w") as fout:
    fout.write(";".join(headers) + "\n")
    for reg_id, d in regions.items():
        c_a = d["depths"].get("ADN_CORIELL_NA24385-A", 0.0)
        c_b = d["depths"].get("ADN_CORIELL_NA24385-B", 0.0)
        c_mean = (c_a + c_b) / 2.0
        
        brca = d["depths"].get("2608120547_BRCA", 0.0)
        mye = d["depths"].get("2608260444_MYELOIDE", 0.0)
        opa_cf = d["depths"].get("2608280855_OPA-CF", 0.0)
        opa = d["depths"].get("2608311668_OPA", 0.0)
        
        max_base = max([brca, mye, opa_cf, opa]) if any([brca, mye, opa_cf, opa]) else 0.0
        ratio = (c_mean / max_base) if max_base > 0 else 0.0
        
        if c_mean == 0:
            diag = "TROU_CAPTURE_TOTAL (0X WES)"
        elif c_mean < 30:
            diag = "COUVERTURE_FAIBLE (<30X)"
        elif c_mean < 100:
            diag = "COUVERTURE_WES_STANDARD (30-100X)"
        else:
            diag = "EXCELLENTE_COUVERTURE (>=100X)"
            
        row = [
            d["gene"], reg_id, d["chr_hg38"], d["start_hg38"], d["end_hg38"], d["status"],
            f"{brca:.1f}", f"{mye:.1f}", f"{opa_cf:.1f}", f"{opa:.1f}",
            f"{max_base:.1f}",
            f"{c_a:.1f}", f"{c_b:.1f}", f"{c_mean:.1f}",
            f"{ratio:.2f}", diag
        ]
        fout.write(";".join(row) + "\n")

print(f"=== Fichier généré avec succès : {csv_file} ===")
EOF
