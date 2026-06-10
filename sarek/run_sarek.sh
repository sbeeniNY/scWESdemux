#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Part 1: Run nf-core/sarek with 3 variant callers (DeepVariant, FreeBayes,
# HaplotypeCaller).
#
# Prerequisites (one-time):
#   # load Nextflow + Singularity however your site provides them, e.g.
#   #   ml nextflow singularity
#   nextflow pull nf-core/sarek -r 3.5.1
#
# Usage:
#   cd sarek
#   bash run_sarek.sh          # full run
#   bash run_sarek.sh -resume  # resume after partial run
#
# Edit the SINGULARITY_* paths below and sarek/nextflow.config before running.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXTRA_ARGS=("$@")

# Load tools (adapt to your environment / module system)
# ml nextflow singularity bcftools htslib

# Redirect Singularity cache/tmp off the home directory to avoid quota issues
export SINGULARITY_CACHEDIR="/path/to/containers/singularity_cache"
export SINGULARITY_TMPDIR="/path/to/tmp/singularity"
export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR}"
mkdir -p "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

cd "${SCRIPT_DIR}"

nextflow run nf-core/sarek \
    -r 3.5.1 \
    -profile singularity \
    -c nextflow.config \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "run_sarek_$(date +%Y%m%d_%H%M%S).log"

bash "${SCRIPT_DIR}/normalize_sarek_vcf_samples.sh"

echo ""
echo "=== Sarek complete ==="
echo "Per-sample VCFs for each caller are in:"
echo "  $(grep -oP 'outdir\s*=\s*"\K[^"]+' nextflow.config)/variant_calling/{deepvariant,freebayes,haplotypecaller}/"
echo ""
echo "Next step: cd .. && snakemake --profile cluster/ --jobs 20"
