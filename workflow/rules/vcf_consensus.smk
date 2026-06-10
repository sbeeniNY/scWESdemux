"""
Multi-caller consensus: for each WES sample, find variant sites called by
>= min_caller_support callers and extract genotypes from the primary caller.
"""

import os

OUTDIR = config["output_dir"]
SAREK_OUT = config["sarek_outdir"]
WES_SAMPLES = list(config["wes_samples"].keys())
CALLERS = config["callers"]
PRIMARY = config["primary_caller"]
MIN_SUPPORT = config["min_caller_support"]


def _sarek_vcf(sample, caller):
    return os.path.join(
        SAREK_OUT,
        "variant_calling",
        caller,
        sample,
        f"{sample}.{caller}.vcf.gz",
    )


rule per_sample_consensus:
    """
    For one WES sample, find sites called by >= min_caller_support callers
    and keep genotypes from the primary caller at those sites.
    """
    input:
        caller_vcfs=lambda wc: [_sarek_vcf(wc.sample, c) for c in CALLERS],
        primary_vcf=lambda wc: _sarek_vcf(wc.sample, PRIMARY),
    output:
        vcf=os.path.join(OUTDIR, "vcf/per_sample/{sample}.consensus.vcf.gz"),
        tbi=os.path.join(OUTDIR, "vcf/per_sample/{sample}.consensus.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/per_sample_consensus/{sample}.log"),
    params:
        min_support=MIN_SUPPORT,
        callers=CALLERS,
    threads: 2
    shell:
        r"""
        set -euo pipefail
        ml bcftools/1.22 htslib/1.21
        mkdir -p $(dirname {output.vcf}) $(dirname {log})

        WORKDIR=$(dirname {output.vcf})/{wildcards.sample}_isec_tmp
        rm -rf "$WORKDIR"
        mkdir -p "$WORKDIR"

        echo "=== Per-sample consensus: {wildcards.sample} ===" > {log}
        echo "Callers: {params.callers}" >> {log}
        echo "Min support: {params.min_support}" >> {log}

        # Extract site positions (CHROM\tPOS) from each caller
        for vcf in {input.caller_vcfs}; do
            caller=$(basename $(dirname $(dirname "$vcf")))
            n=$(bcftools view -H "$vcf" 2>/dev/null | wc -l || echo 0)
            echo "  $caller: $n variants" >> {log}
            bcftools query -f '%CHROM\t%POS\n' "$vcf" 2>/dev/null >> "$WORKDIR/all_sites.txt" || true
        done

        total_sites=$(wc -l < "$WORKDIR/all_sites.txt" || echo 0)
        echo "Total site entries across callers: $total_sites" >> {log}

        # Find sites present in >= min_support callers
        sort "$WORKDIR/all_sites.txt" \
          | uniq -c \
          | awk -v min={params.min_support} '$1 >= min {{ printf "%s\t%s\n", $2, $3 }}' \
          > "$WORKDIR/consensus_sites.tsv"

        n_consensus=$(wc -l < "$WORKDIR/consensus_sites.tsv")
        echo "Consensus sites (>= {params.min_support} callers): $n_consensus" >> {log}

        if [[ "$n_consensus" -eq 0 ]]; then
            echo "WARNING: No consensus sites found. Copying primary caller VCF as-is." >> {log}
            cp {input.primary_vcf} {output.vcf}
        else
            # Filter primary caller VCF to consensus sites only
            bcftools view \
                -T "$WORKDIR/consensus_sites.tsv" \
                -Oz -o {output.vcf} \
                {input.primary_vcf} \
                2>> {log}
        fi

        bcftools index --tbi --force {output.vcf} >> {log} 2>&1

        n_out=$(bcftools view -H {output.vcf} | wc -l)
        echo "Output variants: $n_out" >> {log}

        rm -rf "$WORKDIR"
        """
