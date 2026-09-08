#!/usr/bin/env bash
set -euo pipefail

CACHE_DIR="/NFS/cluster-share/home/mcoquerelle/Explorations/singularity_cache"
mkdir -p "${CACHE_DIR}"

declare -A IMAGES=(
    ["falco_1.2.5.sif"]="bishop/falco:1.2.5"
    ["fastp_0.23.4.sif"]="bishop/fastp:0.23.4"
    ["bwamem2_2.2.1.sif"]="bishop/bwamem2:2.2.1"
    ["bcftools_1.18.sif"]="bishop/bcftools:1.18"
    ["manta_1.6.0.sif"]="bishop/manta:1.6.0"
    ["cnvkit_0.9.10.sif"]="bishop/cnvkit:0.9.10"
    ["gatk_v4.4.0.sif"]="bishop/gatk:v4.4.0"
    ["deepvariant_1.6.0.sif"]="bishop/deepvariant:1.6.0"
    ["vep_110.1.sif"]="bishop/vep:110.1"
    ["mosdepth_0.3.5.sif"]="bishop/mosdepth:0.3.5"
    ["multiqc_1.21.sif"]="bishop/multiqc:1.21"
)

echo "=== Début de la conversion des images SIF pour SLURM ==="

for sif in "${!IMAGES[@]}"; do
    docker_tag="${IMAGES[$sif]}"
    target="${CACHE_DIR}/${sif}"
    tmp="/var/tmp/${sif}"

    if [ -s "${target}" ]; then
        echo ">>> [SKIP] ${sif} existe déjà."
    else
        echo ">>> [BUILD] Création de ${sif} depuis ${docker_tag}..."
        apptainer build --mksquashfs-args "-processors 1" -F "${tmp}" "docker-daemon://${docker_tag}"
        mv "${tmp}" "${target}"
        echo ">>> [SUCCÈS] ${sif} est en place."
    fi
done

echo "=== Toutes les images sont converties avec succès ! ==="
ls -lh "${CACHE_DIR}"/*.sif
