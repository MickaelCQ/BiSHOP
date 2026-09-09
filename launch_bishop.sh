#!/usr/bin/env bash
# ==============================================================================
# Pipeline Exome Constitutionnel BiSHOP V2 — CHU de Nîmes
# Lanceur de production SLURM / Nextflow avec traçabilité et rapports
# ==============================================================================

set -euo pipefail

# Couleurs pour le terminal
RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
BLUE=$'\033[0;34m'
YELLOW=$'\033[1;33m'
BOLD=$'\033[1m'
NC=$'\033[0m'

# Valeurs par défaut
PROFILE="slurm"
INPUT_CSV=""
OUTDIR=""
RESUME=""
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Fonction d'affichage de l'aide
show_help() {
    cat << HELP

${BOLD}==============================================================================${NC}
${BLUE}${BOLD}   BiSHOP V2 — Pipeline Exome Constitutionnel Clinique (CHU de Nîmes)${NC}
${BOLD}==============================================================================${NC}

${BOLD}UTILISATION :${NC}
  ./launch_bishop.sh -i <samplesheet.csv> -o <dossier_sortie> [OPTIONS]

${BOLD}OPTIONS OBLIGATOIRES :${NC}
  ${GREEN}-i, --input${NC} <PATH>       Chemin vers le fichier samplesheet CSV des patients.
  ${GREEN}-o, --outdir${NC} <PATH>      Dossier de destination pour les résultats d'analyse.

${BOLD}OPTIONS FACULTATIVES :${NC}
  ${YELLOW}-p, --profile${NC} <STR>     Profil d'exécution Nextflow (${BOLD}slurm${NC} [défaut] ou ${BOLD}docker${NC}).
  ${YELLOW}-r, --resume${NC}            Active le cache (-resume) pour reprendre un calcul interrompu.
  ${YELLOW}-h, --help${NC}              Affiche cette documentation d'aide.

${BOLD}EXEMPLES D'UTILISATION :${NC}
  ${BLUE}# 1. Exécution de test sur l'échantillon toy en SLURM :${NC}
  ./launch_bishop.sh -i assets/samplesheet_toy_fastq.csv -o results_toy_slurm

  ${BLUE}# 2. Exécution d'une cohorte clinique avec reprise sur incident :${NC}
  ./launch_bishop.sh -i assets/SS_NS2000_All_Fastq_Aout26.csv -o /NFS/cluster-share/projects/WES/Run_2026_09 -r

HELP
    exit 0
}

# Analyse des arguments de la ligne de commande
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_CSV="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME="-resume"
            shift 1
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo -e "${RED}Erreur : Option inconnue '$1'${NC}"
            echo -e "Utilisez ${BOLD}./launch_bishop.sh --help${NC} pour voir la syntaxe."
            exit 1
            ;;
    esac
done

# Vérification des arguments obligatoires
if [[ -z "${INPUT_CSV}" || -z "${OUTDIR}" ]]; then
    echo -e "${RED}Erreur : Les options -i/--input et -o/--outdir sont strictement obligatoires.${NC}"
    echo -e "Tapez ${BOLD}./launch_bishop.sh --help${NC} pour consulter l'aide."
    exit 1
fi

# Vérification de l'existence du fichier CSV
if [[ ! -f "${INPUT_CSV}" ]]; then
    echo -e "${RED}Erreur critique : Le fichier samplesheet introuvable : ${INPUT_CSV}${NC}"
    exit 1
fi

# Création du sous-dossier dédié à la traçabilité Nextflow
chmod -R u+rwX "${OUTDIR}" 2>/dev/null || true
INFO_DIR="${OUTDIR}/pipeline_info"
mkdir -p "${INFO_DIR}"

# Chargement de l'environnement bashrc si nécessaire
source ~/.bashrc 2>/dev/null || true

# Bannière récapitulative
echo -e "${BLUE}${BOLD}==============================================================================${NC}"
echo -e "${BLUE}${BOLD}   LANCEMENT DU PIPELINE BISHOP V2 — ${TIMESTAMP}${NC}"
echo -e "${BLUE}${BOLD}==============================================================================${NC}"
echo -e "  ${BOLD}Samplesheet${NC} : ${INPUT_CSV}"
echo -e "  ${BOLD}Résultats${NC}   : ${OUTDIR}"
echo -e "  ${BOLD}Profil${NC}      : ${PROFILE}"
echo -e "  ${BOLD}Reprise${NC}     : ${RESUME:-Désactivée (From Scratch)}"
echo -e "  ${BOLD}Rapports${NC}    : ${INFO_DIR}/"
echo -e "${BLUE}${BOLD}==============================================================================${NC}"

# Commande Nextflow complète avec tous les rapports d'audit
nextflow run main.nf \
    -profile "${PROFILE}" \
    --input "${INPUT_CSV}" \
    --outdir "${OUTDIR}" \
    ${RESUME} \
    -with-report "${INFO_DIR}/execution_report_${TIMESTAMP}.html" \
    -with-timeline "${INFO_DIR}/timeline_${TIMESTAMP}.html" \
    -with-trace "${INFO_DIR}/trace_${TIMESTAMP}.txt" \
    -with-dag "${INFO_DIR}/dag_${TIMESTAMP}.html"

STATUS=$?

echo -e "${BLUE}${BOLD}==============================================================================${NC}"
if [ ${STATUS} -eq 0 ]; then
    echo -e "${GREEN}${BOLD}   RUN BISHOP V2 TERMINÉ AVEC SUCCÈS !${NC}"
    echo -e "   Consultez le rapport d'exécution : ${INFO_DIR}/execution_report_${TIMESTAMP}.html"
else
    echo -e "${RED}${BOLD}   ÉCHEC DU RUN BISHOP V2 (Code de sortie : ${STATUS})${NC}"
fi
echo -e "${BLUE}${BOLD}==============================================================================${NC}"

exit ${STATUS}
