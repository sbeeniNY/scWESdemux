#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Caller overlap summary: 3-way Venn counts per WES sample
#
# Usage:
#   bash caller_overlap_summary.sh <sarek_outdir> [output_tsv]
#
# Example:
#   bash caller_overlap_summary.sh /path/to/output/sarek_results
#
# Samples are auto-detected from the haplotypecaller output subdirectories.
# =============================================================================

SAREK_OUT="${1:?Usage: $0 <sarek_outdir> [output_tsv]}"
OUTPUT="${2:-${SAREK_OUT}/caller_overlap_summary.tsv}"

ml bcftools htslib 2>/dev/null || true

# Auto-detect sample names from one caller's output directory.
HC_DIR="${SAREK_OUT}/variant_calling/haplotypecaller"
if [[ ! -d "${HC_DIR}" ]]; then
    echo "ERROR: ${HC_DIR} not found. Has sarek finished?" >&2
    exit 1
fi
mapfile -t SAMPLES < <(find "${HC_DIR}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: No sample directories under ${HC_DIR}" >&2
    exit 1
fi

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

header="sample\tDV_total\tFB_total\tHC_total\tDV_only\tFB_only\tHC_only\tDV_FB\tDV_HC\tFB_HC\tall_three\tconsensus_ge2"
echo -e "$header" > "$OUTPUT"

for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample ..."

    DV="${SAREK_OUT}/variant_calling/deepvariant/${sample}/${sample}.deepvariant.vcf.gz"
    FB="${SAREK_OUT}/variant_calling/freebayes/${sample}/${sample}.freebayes.vcf.gz"
    HC="${SAREK_OUT}/variant_calling/haplotypecaller/${sample}/${sample}.haplotypecaller.vcf.gz"

    for f in "$DV" "$FB" "$HC"; do
        if [[ ! -f "$f" ]]; then
            echo "  WARNING: $f not found, skipping $sample"
            continue 2
        fi
    done

    # Extract CHROM:POS from each caller
    bcftools query -f '%CHROM:%POS\n' "$DV" 2>/dev/null | sort -u > "$TMPDIR/${sample}_dv.txt"
    bcftools query -f '%CHROM:%POS\n' "$FB" 2>/dev/null | sort -u > "$TMPDIR/${sample}_fb.txt"
    bcftools query -f '%CHROM:%POS\n' "$HC" 2>/dev/null | sort -u > "$TMPDIR/${sample}_hc.txt"

    dv_total=$(wc -l < "$TMPDIR/${sample}_dv.txt")
    fb_total=$(wc -l < "$TMPDIR/${sample}_fb.txt")
    hc_total=$(wc -l < "$TMPDIR/${sample}_hc.txt")

    # Pairwise intersections
    dv_fb_sites=$(comm -12 "$TMPDIR/${sample}_dv.txt" "$TMPDIR/${sample}_fb.txt" | wc -l)
    dv_hc_sites=$(comm -12 "$TMPDIR/${sample}_dv.txt" "$TMPDIR/${sample}_hc.txt" | wc -l)
    fb_hc_sites=$(comm -12 "$TMPDIR/${sample}_fb.txt" "$TMPDIR/${sample}_hc.txt" | wc -l)

    # 3-way intersection
    comm -12 "$TMPDIR/${sample}_dv.txt" "$TMPDIR/${sample}_fb.txt" > "$TMPDIR/${sample}_dv_fb.txt"
    all3=$(comm -12 "$TMPDIR/${sample}_dv_fb.txt" "$TMPDIR/${sample}_hc.txt" | wc -l)

    # Exclusive pairwise (subtract 3-way)
    dv_fb=$((dv_fb_sites - all3))
    dv_hc=$((dv_hc_sites - all3))
    fb_hc=$((fb_hc_sites - all3))

    # Exclusive singles
    dv_only=$((dv_total - dv_fb - dv_hc - all3))
    fb_only=$((fb_total - dv_fb - fb_hc - all3))
    hc_only=$((hc_total - dv_hc - fb_hc - all3))

    # Consensus: sites in >= 2 callers
    consensus=$((dv_fb + dv_hc + fb_hc + all3))

    echo -e "${sample}\t${dv_total}\t${fb_total}\t${hc_total}\t${dv_only}\t${fb_only}\t${hc_only}\t${dv_fb}\t${dv_hc}\t${fb_hc}\t${all3}\t${consensus}" >> "$OUTPUT"
done

echo ""
echo "=== Caller Overlap Summary ==="
column -t "$OUTPUT"
echo ""
echo "Written to: $OUTPUT"
echo ""
echo "Legend:"
echo "  DV = DeepVariant, FB = FreeBayes, HC = HaplotypeCaller"
echo "  DV_FB = shared by DV+FB only (not HC), etc."
echo "  all_three = called by all 3 callers"
echo "  consensus_ge2 = sites called by >= 2 callers (used in pipeline)"
