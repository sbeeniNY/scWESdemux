#!/usr/bin/env bash
set -euo pipefail
# Normalize nf-core/sarek VCF sample names from <patient>_<sample>
# back to the logical donor/sample ID used by the demux config.
# Handles the 3 callers: deepvariant, freebayes, haplotypecaller.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLESHEET="${SCRIPT_DIR}/samplesheet.csv"
CONFIG="${SCRIPT_DIR}/nextflow.config"

CALLERS=("deepvariant" "freebayes" "haplotypecaller")

if command -v ml >/dev/null 2>&1; then
    ml bcftools htslib 2>/dev/null || true
fi

OUTDIR="$(grep -oP 'outdir\s*=\s*"\K[^"]+' "${CONFIG}")"
if [[ -z "${OUTDIR}" ]]; then
    echo "ERROR: Could not parse Sarek outdir from ${CONFIG}" >&2
    exit 1
fi

if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "ERROR: Missing samplesheet: ${SAMPLESHEET}" >&2
    exit 1
fi

echo "Normalizing VCF sample names in: ${OUTDIR}"

for caller in "${CALLERS[@]}"; do
    echo ""
    echo "--- Caller: ${caller} ---"

    tail -n +2 "${SAMPLESHEET}" | while IFS=, read -r patient sample lane fastq_1 fastq_2; do
        [[ -n "${sample}" ]] || continue

        old_name="${patient}_${sample}"
        new_name="${sample}"
        vcf="${OUTDIR}/variant_calling/${caller}/${sample}/${sample}.${caller}.vcf.gz"

        if [[ ! -f "${vcf}" ]]; then
            echo "  WARN: VCF not found (skipping): ${vcf}"
            continue
        fi

        current_samples="$(bcftools query -l "${vcf}")"

        if grep -Fxq "${new_name}" <<< "${current_samples}" && ! grep -Fxq "${old_name}" <<< "${current_samples}"; then
            echo "  ${caller}/${sample}: already normalized (${new_name})"
            bcftools index --tbi --force "${vcf}"
            continue
        fi

        if ! grep -Fxq "${old_name}" <<< "${current_samples}"; then
            echo "ERROR: ${vcf} does not contain expected sample ${old_name}" >&2
            echo "Current samples:" >&2
            printf '%s\n' "${current_samples}" >&2
            exit 1
        fi

        work_dir="$(dirname "${vcf}")"
        rename_file="${work_dir}/sample_rename.tsv"
        tmp_vcf="${work_dir}/${sample}.${caller}.renamed.tmp.vcf.gz"

        printf '%s\t%s\n' "${old_name}" "${new_name}" > "${rename_file}"
        bcftools reheader -s "${rename_file}" -o "${tmp_vcf}" "${vcf}"
        bcftools index --tbi --force "${tmp_vcf}"

        renamed_samples="$(bcftools query -l "${tmp_vcf}")"
        if ! grep -Fxq "${new_name}" <<< "${renamed_samples}" || grep -Fxq "${old_name}" <<< "${renamed_samples}"; then
            echo "ERROR: Failed to normalize sample name in ${tmp_vcf}" >&2
            exit 1
        fi

        mv -f "${tmp_vcf}" "${vcf}"
        mv -f "${tmp_vcf}.tbi" "${vcf}.tbi"
        rm -f "${rename_file}"

        echo "  ${caller}/${sample}: ${old_name} -> ${new_name}"
    done
done

echo ""
echo "VCF sample-name normalization complete for all callers."
