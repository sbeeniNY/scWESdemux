"""
Genetic demultiplexing per pool using Souporcell.

Two mutually exclusive modes (Souporcell does NOT allow both simultaneously):

1. common_variants mode (recommended):
   --common_variants only (de novo clustering using genome-wide common SNPs).
   After Souporcell finishes, assign_clusters.py matches numeric cluster IDs
   to donor names via genotype concordance against the WES donor VCF.

2. known_genotypes mode (WES-only fallback):
   --known_genotypes + --known_genotypes_sample_names (supervised, WES sites only).
   Souporcell assigns donor names directly.

Mode is selected by config['souporcell']['common_variants']:
  - Non-empty path -> common_variants mode
  - Empty string   -> known_genotypes mode
"""

import os

OUTDIR = config["output_dir"]
POOLS = list(config["pools"].keys())
FILTER_ENABLED = config["vcf_filter"]["enabled"]
INFORMATIVE_FILTER = config.get("informative_filter", {}).get("enabled", True)
COMMON_VARIANTS = config["souporcell"].get("common_variants", "")

if FILTER_ENABLED:
    _COHORT_VCF = os.path.join(OUTDIR, "vcf/cohort.pass.vcf.gz")
    _COHORT_TBI = os.path.join(OUTDIR, "vcf/cohort.pass.vcf.gz.tbi")
else:
    _COHORT_VCF = os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz")
    _COHORT_TBI = os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz.tbi")


def _pool_bam(wc):
    cr_dir = config["pools"][wc.pool]["cellranger_dir"]
    cr_id = os.path.basename(os.path.dirname(cr_dir))
    return os.path.join(cr_dir, "per_sample_outs", cr_id, "sample_alignments.bam")


def _pool_barcodes(wc):
    return os.path.join(
        config["pools"][wc.pool]["cellranger_dir"],
        "filtered_feature_bc_matrix",
        "barcodes.tsv.gz",
    )


def _pool_donors_csv(pool):
    return ",".join(config["pools"][pool]["donors"])


rule subset_pool_donor_vcf:
    input:
        vcf=_COHORT_VCF,
        tbi=_COHORT_TBI,
    output:
        vcf=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz"),
        tbi=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/subset_pool_donor_vcf/{pool}.log"),
    threads: 2
    params:
        donors=lambda wc: _pool_donors_csv(wc.pool),
        informative_filter=INFORMATIVE_FILTER,
    shell:
        r"""
        set -euo pipefail
        ml bcftools/1.22 htslib/1.21

        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        : > {log}

        echo "Requested donors: {params.donors}" >> {log}
        echo "Informative filter: {params.informative_filter}" >> {log}
        echo "Available VCF samples:" >> {log}
        bcftools query -l {input.vcf} >> {log}

        if ! bcftools view -h {input.vcf} | grep -q '^##FORMAT=<ID=GT,'; then
            echo "ERROR: input VCF has no FORMAT/GT header; cannot subset donor genotypes." >> {log}
            exit 1
        fi

        missing=0
        for donor in $(echo {params.donors} | tr ',' ' '); do
            if ! bcftools query -l {input.vcf} | grep -Fxq "$donor"; then
                echo "ERROR: requested donor not found in VCF header: $donor" >> {log}
                missing=1
            fi
        done
        if [[ "$missing" -ne 0 ]]; then
            exit 1
        fi

        if [[ "{params.informative_filter}" == "True" ]]; then
            # Subset donors + remove missing GT + keep only donor-discriminating sites
            bcftools view -s {params.donors} -a {input.vcf} 2>> {log} \
              | bcftools view -e 'GT[*]="mis"' 2>> {log} \
              | bcftools view -e 'COUNT(GT="RR")=N_SAMPLES || COUNT(GT="het")=N_SAMPLES || COUNT(GT="AA")=N_SAMPLES' \
                  -Oz -o {output.vcf} 2>> {log}
        else
            # Subset donors + remove missing GT only (no informative filter)
            bcftools view -s {params.donors} -a {input.vcf} 2>> {log} \
              | bcftools view -e 'GT[*]="mis"' -Oz -o {output.vcf} 2>> {log}
        fi

        bcftools index --tbi --force {output.vcf} >> {log} 2>&1

        n_sites=$(bcftools view -H {output.vcf} | wc -l)
        echo "Final variant sites for Souporcell: $n_sites" >> {log}
        """


rule souporcell:
    input:
        bam=lambda wc: _pool_bam(wc),
        barcodes=lambda wc: _pool_barcodes(wc),
        donor_vcf=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz"),
        donor_tbi=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        outdir=directory(os.path.join(OUTDIR, "demux/{pool}/souporcell")),
        clusters=os.path.join(OUTDIR, "demux/{pool}/souporcell/clusters.tsv"),
    log:
        os.path.join(OUTDIR, "logs/souporcell/{pool}.log"),
    threads: config["souporcell"]["threads"]
    params:
        sif=config["souporcell"]["sif"],
        n_donor=lambda wc: config["pools"][wc.pool]["n_donor"],
        donor_names=lambda wc: " ".join(config["pools"][wc.pool]["donors"]),
        common_variants=COMMON_VARIANTS,
        extra=config["souporcell"]["extra"],
        assign_script=os.path.join(os.path.dirname(workflow.snakefile), "scripts", "assign_clusters.py"),
    shell:
        r"""
        set -euo pipefail
        ml singularity/3.11.0 htslib/1.21
        mkdir -p $(dirname {log})
        mkdir -p {output.outdir}
        mkdir -p {output.outdir}/logs

        : > {log}

        decompressed_vcf={output.outdir}/donor_genotype.vcf
        bgzip -dkc {input.donor_vcf} > "$decompressed_vcf"

        if [[ -n "{params.common_variants}" ]]; then
            echo "MODE: common_variants (de novo + post-hoc donor matching)" >> {log}
            echo "Common variants panel: {params.common_variants}" >> {log}

            singularity exec \
                --bind /sc:/sc \
                {params.sif} \
                souporcell_pipeline.py \
                    -i {input.bam} \
                    -b {input.barcodes} \
                    -f {input.ref} \
                    -t {threads} \
                    -o {output.outdir} \
                    -k {params.n_donor} \
                    --common_variants {params.common_variants} \
                    {params.extra} \
                    >> {log} 2>&1

            echo "Post-hoc donor assignment via WES genotype concordance..." >> {log}
            python3 {params.assign_script} \
                --cluster-vcf {output.outdir}/cluster_genotypes.vcf \
                --donor-vcf "$decompressed_vcf" \
                --clusters-tsv {output.outdir}/clusters.tsv \
                --report {output.outdir}/donor_assignment_report.tsv \
                >> {log} 2>&1
        else
            echo "MODE: known_genotypes (supervised, WES-only)" >> {log}

            singularity exec \
                --bind /sc:/sc \
                {params.sif} \
                souporcell_pipeline.py \
                    -i {input.bam} \
                    -b {input.barcodes} \
                    -f {input.ref} \
                    -t {threads} \
                    -o {output.outdir} \
                    -k {params.n_donor} \
                    --known_genotypes "$decompressed_vcf" \
                    --known_genotypes_sample_names {params.donor_names} \
                    {params.extra} \
                    >> {log} 2>&1
        fi
        """


rule demux_summary:
    input:
        clusters=expand(
            os.path.join(OUTDIR, "demux/{pool}/souporcell/clusters.tsv"),
            pool=POOLS,
        ),
        config_yaml="config/config.yaml",
    output:
        tsv=os.path.join(OUTDIR, "demux/demux_summary.tsv"),
    log:
        os.path.join(OUTDIR, "logs/demux_summary/all.log"),
    threads: 1
    script:
        "../scripts/demux_summary.py"
