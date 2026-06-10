"""
VCF merge: combine per-sample consensus VCFs into a multi-sample cohort VCF.
Optional quality filter controlled by config['vcf_filter'].
"""

import os

OUTDIR = config["output_dir"]
WES_SAMPLES = list(config["wes_samples"].keys())
FILTER_ENABLED = config["vcf_filter"]["enabled"]
QUAL_MIN = config["vcf_filter"]["qual_min"]


rule merge_consensus_vcfs:
    input:
        vcfs=expand(
            os.path.join(OUTDIR, "vcf/per_sample/{sample}.consensus.vcf.gz"),
            sample=WES_SAMPLES,
        ),
    output:
        vcf=os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz"),
        tbi=os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/merge_consensus_vcfs/cohort.log"),
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        ml bcftools/1.22 htslib/1.21
        bcftools merge \
            --merge both \
            --threads {threads} \
            -Oz -o {output.vcf} \
            {input.vcfs} \
            > {log} 2>&1
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """


rule filter_pass:
    input:
        vcf=os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz"),
        tbi=os.path.join(OUTDIR, "vcf/cohort.consensus.vcf.gz.tbi"),
    output:
        vcf=os.path.join(OUTDIR, "vcf/cohort.pass.vcf.gz"),
        tbi=os.path.join(OUTDIR, "vcf/cohort.pass.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/filter_pass/cohort.log"),
    params:
        qual_min=QUAL_MIN,
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) $(dirname {log})
        ml bcftools/1.22 htslib/1.21
        bcftools view \
            -f 'PASS,.' \
            -i 'QUAL>={params.qual_min}' \
            -Oz -o {output.vcf} \
            {input.vcf} \
            2> {log}
        bcftools index --tbi --force {output.vcf} >> {log} 2>&1
        """
