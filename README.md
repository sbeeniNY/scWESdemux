# scWESdemux_SarekSoup
## WES -> nf-core/sarek (Multi-caller Consensus) -> Souporcell (v4)

**WES-guided demultiplexing for pooled single-cell RNA-seq using nf-core/sarek, multi-caller consensus variants, and Souporcell.**

`scWESdemux` is a two-part workflow for assigning pooled scRNA-seq barcodes to donors when matched whole-exome sequencing (WES) data are available. It first calls donor variants from WES using `nf-core/sarek`, builds a multi-caller consensus VCF, and then runs `Souporcell` on pooled scRNA-seq BAMs. In the default mode, Souporcell uses a genome-wide common SNP panel for de novo clustering, and `scWESdemux` assigns Souporcell cluster IDs back to donor identities by genotype concordance against the WES-derived donor VCF.

## Unique features: Why scWESdemux pipeline? 

WES-only demultiplexing is limited because exome variants are concentrated in coding regions, whereas 10x 3' scRNA-seq reads are biased toward 3' UTRs. In practice, only a small fraction of WES variants are covered by scRNA-seq reads. `scWESdemux` addresses this by combining:

- **WES-derived donor genotypes** for donor identity matching
- **Multi-caller WES consensus variants** from DeepVariant, FreeBayes, and GATK HaplotypeCaller
- **Souporcell `--common_variants` mode** using genome-wide common SNPs for de novo clustering from scRNA-seq reads
- **Post-hoc cluster-to-donor assignment** by genotype concordance
- **Configurable VCF and informative-site filters** for benchmarking and sensitivity analysis

This design is useful when matched WES data exist but WES-only variant sites provide insufficient coverage for confident scRNA-seq demultiplexing.

## Workflow overview

```text
Matched WES FASTQs
        │
        ▼
Part 1: nf-core/sarek
  - BWA-MEM2 alignment
  - duplicate marking
  - BQSR
  - variant calling with DeepVariant, FreeBayes, HaplotypeCaller
        │
        ▼
Per-donor VCFs from three callers
        │
        ▼
Part 2: Snakemake workflow
  - normalize donor sample names
  - build per-donor multi-caller consensus VCFs
  - merge donor VCFs into cohort VCF
  - subset donor VCFs per scRNA-seq pool
  - optionally retain donor-discriminating sites only
  - run Souporcell with common variants and skip_remap
  - assign numeric clusters to donor IDs
  - summarize singlet, doublet, and unassigned calls
```

## Pipeline structure

```text
scWESdemux/
├── config/
│   └── config.yaml                     # Central configuration file
├── sarek/
│   ├── samplesheet.csv                 # nf-core/sarek input samplesheet
│   ├── nextflow.config                 # Sarek configuration
│   ├── run_sarek.sh                    # Launch script for WES variant calling
│   └── normalize_sarek_vcf_samples.sh  # Rename Sarek VCF samples to donor IDs
├── workflow/
│   ├── Snakefile                       # Snakemake entry point
│   ├── rules/
│   │   ├── vcf_consensus.smk           # Multi-caller consensus per donor
│   │   ├── vcf_merge.smk               # Cohort merge and filtering
│   │   └── demux.smk                   # Souporcell and summary rules
│   └── scripts/
│       ├── assign_clusters.py          # Cluster-to-donor assignment
│       └── demux_summary.py            # Aggregate demultiplexing summary
├── scripts/
│   ├── caller_overlap_summary.sh       # Caller overlap QC
│   └── demux_summary_standalone.py     # Standalone summary utility
├── cluster/
│   ├── config.yaml                     # Snakemake cluster profile example
│   ├── lsf_submit.sh                   # LSF submission wrapper
│   └── lsf_status.sh                   # LSF job status helper
└── README.md
```

--- 

## Requirements

### Software

- Nextflow
- nf-core/sarek
- Snakemake
- Singularity or Apptainer
- bcftools
- htslib/tabix
- Python 3
- Souporcell container

The workflow was developed for an HPC environment using LSF, but the Snakemake profile can be adapted to other schedulers.

### Reference files

Required reference inputs are configured in `config/config.yaml`.

1. **Reference genome FASTA**  
   Use the same reference build and chromosome naming convention for WES and scRNA-seq. For Cell Ranger BAMs, this usually means using the FASTA from the matching Cell Ranger reference directory.

2. **WES target BED file**  
   For WES data, provide the target BED file used by the capture panel. The original implementation used the Twist Bioscience Human Comprehensive Exome covered targets BED for hg38.

3. **Common variants panel**  
   A genome-wide common SNP panel is used by Souporcell through `--common_variants`. The original implementation used the cellSNP-lite 1000 Genomes Phase 3 common SNP list for hg38 with chr-prefixed chromosomes.

4. **Souporcell container**  
   Provide a Singularity/Apptainer image path in `config/config.yaml`.

---

## Input data

Edit `config/config.yaml` to define:

- WES FASTQ files for each donor
- scRNA-seq Cell Ranger `outs/` directory for each pool
- donor IDs included in each pool
- number of expected donors per pool
- reference genome FASTA
- WES target BED
- Souporcell container
- common variants VCF
- variant caller settings
- VCF filtering options
- output directory

Example pool configuration:

```yaml
pools:
  Pool1:
    cellranger_dir: "/path/to/Pool1/cellranger/outs"
    donors: ["DonorA", "DonorB", "DonorC"]
    n_donor: 3
  Pool2:
    cellranger_dir: "/path/to/Pool2/cellranger/outs"
    donors: ["DonorD", "DonorE"]
    n_donor: 2
```

## Quick start

### 1. Configure WES variant calling

Edit:

```text
sarek/samplesheet.csv
sarek/nextflow.config
config/config.yaml
```

Make sure donor names are consistent across WES FASTQs, pool definitions, and downstream expected donor IDs.

### 2. Run nf-core/sarek

```bash
cd sarek
bash run_sarek.sh
```

To resume an interrupted Sarek run:

```bash
bash run_sarek.sh -resume
```

Sarek generates variant calls for each donor and each caller:

```text
{sarek_outdir}/variant_calling/{deepvariant,freebayes,haplotypecaller}/{sample}/{sample}.{caller}.vcf.gz
```

After Sarek completes, `normalize_sarek_vcf_samples.sh` reheaders VCF sample names from Sarek sample conventions to donor IDs used by the demultiplexing workflow.

### 3. Check caller overlap

Before running demultiplexing, inspect how much the three WES callers agree:

```bash
bash scripts/caller_overlap_summary.sh /path/to/sarek_results
```

This produces a per-sample TSV with caller-specific and consensus-site counts:

```text
sample  DV_total  FB_total  HC_total  DV_only  FB_only  HC_only  DV_FB  DV_HC  FB_HC  all_three  consensus_ge2
```

`consensus_ge2` is the number of sites called by at least two of the three callers.

### 4. Run consensus VCF generation and Souporcell

From the repository root:

```bash
snakemake -n --profile cluster/
snakemake --profile cluster/ --jobs 20
```

To resume:

```bash
snakemake --profile cluster/ --jobs 20 --rerun-incomplete
```

The Snakemake workflow automatically:

1. Creates per-donor consensus VCFs using sites called by at least two callers
2. Uses HaplotypeCaller genotypes at consensus sites
3. Merges per-donor consensus VCFs into a cohort VCF
4. Optionally applies PASS/QUAL filtering
5. Subsets the VCF for each pool using only the donors present in that pool
6. Optionally keeps only donor-discriminating sites
7. Runs Souporcell with `--common_variants` and `--skip_remap TRUE`
8. Assigns Souporcell numeric clusters to donor IDs by WES genotype concordance
9. Aggregates demultiplexing metrics into `demux_summary.tsv`

## Key configuration options

### Variant callers

```yaml
callers:
  - deepvariant
  - freebayes
  - haplotypecaller

primary_caller: haplotypecaller
min_caller_support: 2
```

The default consensus strategy keeps sites called by at least two of the three callers and uses HaplotypeCaller genotypes for retained sites.

### VCF quality filtering

```yaml
vcf_filter:
  enabled: true
  qual_min: 20
```

When enabled, Souporcell receives a PASS/QUAL-filtered consensus VCF. When disabled, Souporcell receives the raw multi-caller consensus VCF.

### Donor-discriminating site filtering

```yaml
informative_filter:
  enabled: true
```

When enabled, the per-pool donor VCF keeps only sites where at least one donor differs from another donor in the same pool. Sites where all donors share the same genotype are removed because they do not help distinguish donor identities.

### Souporcell common variants

```yaml
souporcell:
  sif: "/path/to/containers/souporcell_latest.sif"
  threads: 8
  common_variants: "/path/to/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
  extra: "--skip_remap TRUE"   # passed verbatim to souporcell_pipeline.py
```

`--common_variants` and `--known_genotypes` are mutually exclusive in Souporcell. `scWESdemux` therefore uses Souporcell common-variants mode for de novo clustering, then maps clusters to donor IDs after the run.

## Outputs

```
{output_dir}/
├── vcf/
│   ├── per_sample/
│   │   └── {sample}.consensus.vcf.gz       # per-sample multi-caller consensus
│   ├── cohort.consensus.vcf.gz             # merged consensus (all donors)
│   └── cohort.pass.vcf.gz                  # QUAL>=20 filtered (if enabled)
├── demux/
│   ├── {pool}/
│   │   ├── donor_genotype.vcf.gz           # per-pool donor subset (+ .tbi)
│   │   └── souporcell/
│   │       ├── clusters.tsv                # barcode -> assignment (final, post-troublet)
│   │       ├── clusters_with_samples.tsv   # clusters.tsv + donor_id + sample_label columns
│   │       ├── cluster_genotypes.vcf       # Souporcell-inferred cluster genotypes
│   │       ├── donor_genotype.vcf          # decompressed per-pool donor VCF
│   │       ├── ambient_rna.txt             # estimated ambient RNA fraction
│   │       ├── donor_assignment_report.tsv # concordance matrix (common_variants mode ONLY)
│   │       ├── clusters_tmp.tsv            # Souporcell internal (pre-troublet)
│   │       ├── common_variants_covered.vcf # Souporcell internal (covered sites)
│   │       ├── alt.mtx / ref.mtx           # vartrix allele-count matrices
│   │       ├── barcodes.tsv                # barcodes used
│   │       ├── depth_merged.bed            # per-site depth
│   │       ├── vartrix.out                 # vartrix log
│   │       ├── *.done                      # stage marker files
│   │       └── logs/                       # per-stage logs
│   └── demux_summary.tsv                   # cohort-wide summary
└── sarek_results/                          # nf-core/sarek output tree
```

**Key per-pool files:**

- **`clusters_with_samples.tsv`** — the human-readable demux result. It is
  `clusters.tsv` plus two added columns (`sample_assignment` = donor ID,
  `sample_label` = final label from `config.yaml > sample_labels`). Written by
  `demux_summary.py` for every pool.
- **`clusters.tsv`** — raw Souporcell output. In **common_variants mode** the
  `assignment` column is rewritten from numeric cluster IDs to donor names by
  `assign_clusters.py`. In **known_genotypes mode** Souporcell writes the
  assignment directly.
- **`donor_assignment_report.tsv`** — cluster-to-donor concordance matrix. Only
  produced in **common_variants mode** (it is the output of `assign_clusters.py`).
  Absent in known_genotypes (WES-only) runs.

The top-level summary file reports donor-level demultiplexing metrics:

```text
pool   donor_id   mapped_label   n_cells   doublet_rate   unassigned_rate   total_cells
```

---

## Benchmarking filtered vs unfiltered VCFs

To compare different VCF filtering strategies, rerun only Part 2 with a new output directory.

Filtered run:

```yaml
vcf_filter:
  enabled: true
  qual_min: 20
output_dir: "/path/to/output_filtered"
```

Unfiltered run:

```yaml
vcf_filter:
  enabled: false
output_dir: "/path/to/output_unfiltered"
```

Then rerun:

```bash
snakemake --profile cluster/ --jobs 20
```

Because Sarek output is shared, only the downstream consensus, filtering, Souporcell, and summary steps need to be rerun.

---

## Citation

If you use `scWESdemux`, please cite this repository and the underlying tools and resources used in your analysis. At minimum, cite the workflow, `nf-core/sarek`, `Souporcell`, the variant callers used for WES variant discovery, and the common variant resource used for demultiplexing.

### scWESdemux

If this repository has a Zenodo DOI or tagged release, cite that release. Until then, cite the GitHub repository:

```bibtex
@software{scwesdemux,
  title        = {scWESdemux: WES-guided demultiplexing for pooled single-cell RNA-seq},
  author       = {Cho, Subin},
  year         = {2026},
  url          = {https://github.com/sbeeniNY/scWESdemux},
  version      = {v.1}
}
```

### Core workflow and demultiplexing

- Hanssen F, Garcia MU, Folkersen L, Pedersen AS, Lescai F, Jodoin S, Miller E, Wacker O, Smith N, nf-core community, Gabernet G, Nahnsen S. **Scalable and efficient DNA sequencing analysis on different compute infrastructures aiding variant discovery.** *NAR Genomics and Bioinformatics*. 2024;6(2):lqae031. doi: [10.1093/nargab/lqae031](https://doi.org/10.1093/nargab/lqae031).  
  Used for: `nf-core/sarek`.

- Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **The nf-core framework for community-curated bioinformatics pipelines.** *Nature Biotechnology*. 2020;38:276–278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).  
  Used for: the nf-core framework.

- Heaton H, Talman AM, Knights A, Imaz M, Gaffney DJ, Durbin R, Hemberg M, Lawniczak MKN. **Souporcell: robust clustering of single-cell RNA-seq data by genotype without reference genotypes.** *Nature Methods*. 2020;17:615–620. doi: [10.1038/s41592-020-0820-1](https://doi.org/10.1038/s41592-020-0820-1).  
  Used for: genotype-based scRNA-seq demultiplexing, doublet detection, and ambient RNA estimation.

### Workflow engines and execution environment

- Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. **Nextflow enables reproducible computational workflows.** *Nature Biotechnology*. 2017;35:316–319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820).  
  Used for: running `nf-core/sarek`.

- Köster J, Rahmann S. **Snakemake—a scalable bioinformatics workflow engine.** *Bioinformatics*. 2012;28(19):2520–2522. doi: [10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480).  
  Used for: downstream consensus VCF generation, per-pool VCF subsetting, Souporcell execution, and summary aggregation.

- Kurtzer GM, Sochat V, Bauer MW. **Singularity: Scientific containers for mobility of compute.** *PLOS ONE*. 2017;12(5):e0177459. doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459).  
  Used for: containerized execution of Souporcell and/or pipeline components on HPC systems.

### Alignment, variant calling, and VCF processing

- Li H. **Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.** *arXiv*. 2013. [arXiv:1303.3997](https://arxiv.org/abs/1303.3997).  
  Used for: BWA-MEM/BWA-MEM2 short-read alignment in `nf-core/sarek`.

- Li H, Durbin R. **Fast and accurate short read alignment with Burrows-Wheeler transform.** *Bioinformatics*. 2009;25(14):1754–1760. doi: [10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324).  
  Used for: the original BWA aligner.

- Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, et al. **A universal SNP and small-indel variant caller using deep neural networks.** *Nature Biotechnology*. 2018;36:983–987. doi: [10.1038/nbt.4235](https://doi.org/10.1038/nbt.4235).  
  Used for: `DeepVariant`.

- Garrison E, Marth G. **Haplotype-based variant detection from short-read sequencing.** *arXiv*. 2012. [arXiv:1207.3907](https://arxiv.org/abs/1207.3907).  
  Used for: `FreeBayes`.

- Poplin R, Ruano-Rubio V, DePristo MA, Fennell TJ, Carneiro MO, Van der Auwera GA, et al. **Scaling accurate genetic variant discovery to tens of thousands of samples.** *bioRxiv*. 2018. doi: [10.1101/201178](https://doi.org/10.1101/201178).  
  Used for: GATK `HaplotypeCaller` and joint germline variant discovery.

- Van der Auwera GA, O'Connor BD. **Genomics in the Cloud: Using Docker, GATK, and WDL in Terra.** 1st ed. O'Reilly Media; 2020.  
  Used for: general GATK best-practice citation.

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. **Twelve years of SAMtools and BCFtools.** *GigaScience*. 2021;10(2):giab008. doi: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).  
  Used for: `bcftools`, `tabix`, and HTSlib-based VCF/BCF processing.

### Common variant resources

- The 1000 Genomes Project Consortium. **A global reference for human genetic variation.** *Nature*. 2015;526:68–74. doi: [10.1038/nature15393](https://doi.org/10.1038/nature15393).  
  Used for: 1000 Genomes common variant resources.

- Huang X, Huang Y. **Cellsnp-lite: an efficient tool for genotyping single cells.** *Bioinformatics*. 2021;37(23):4569–4571. doi: [10.1093/bioinformatics/btab358](https://doi.org/10.1093/bioinformatics/btab358).  
  Used if you use the cellSNP-lite curated 1000 Genomes common SNP list, for example `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz`.
