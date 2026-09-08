# BiSHOP V2 - Pipeline Exome Constitutionnel Multi-Activités

Pipeline DSL2 Nextflow pour la biologie moléculaire hospitalière.

## Exécution
```bash
nextflow run main.nf \
  -profile singularity \
  --input assets/samplesheet.csv \
  --fasta /chemin/GRCh38.fa \
  --fasta_fai /chemin/GRCh38.fa.fai \
  --bwa_index /chemin/bwa-mem2_index/GRCh38 \
  --vep_cache /chemin/vep_cache \
  --outdir ./results
