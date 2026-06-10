#!/usr/bin/env python3
"""
Post-hoc cluster-to-donor assignment using genotype concordance.

When Souporcell runs with --common_variants (de novo mode), clusters are
numbered 0, 1, 2, ... This script matches those cluster IDs to actual
donor names by comparing cluster_genotypes.vcf (Souporcell output) against
a WES-derived donor_genotype.vcf via genotype concordance.

Rewrites clusters.tsv in place with donor names replacing numeric indices.

Usage:
    python assign_clusters.py \
        --cluster-vcf <souporcell_outdir>/cluster_genotypes.vcf \
        --donor-vcf <donor_genotype.vcf> \
        --clusters-tsv <souporcell_outdir>/clusters.tsv \
        --report <souporcell_outdir>/donor_assignment_report.tsv
"""

import argparse
import csv
import os
import sys


def parse_vcf_genotypes(vcf_path):
    """Parse VCF and return (samples, {(chrom, pos): {sample: GT_tuple}})."""
    genotypes = {}
    samples = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                samples = fields[9:]
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = fields[1]
            alt = fields[4]

            if "," in alt:
                continue

            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            site_gts = {}
            for i, sample in enumerate(samples):
                sample_data = fields[9 + i].split(":")
                if gt_idx >= len(sample_data):
                    continue
                gt_str = sample_data[gt_idx]

                sep = "/" if "/" in gt_str else "|"
                alleles = gt_str.split(sep)
                if "." in alleles:
                    continue
                try:
                    gt = tuple(sorted(int(a) for a in alleles))
                except ValueError:
                    continue
                site_gts[sample] = gt

            if site_gts:
                genotypes[(chrom, pos)] = site_gts

    return samples, genotypes


def compute_concordance(cluster_gts, donor_gts):
    """Compute genotype concordance for each cluster-donor pair at shared sites."""
    shared_positions = set(cluster_gts.keys()) & set(donor_gts.keys())

    cluster_samples = set()
    donor_samples = set()
    for pos in shared_positions:
        cluster_samples.update(cluster_gts[pos].keys())
        donor_samples.update(donor_gts[pos].keys())

    cluster_samples = sorted(cluster_samples, key=lambda x: int(x) if x.isdigit() else x)
    donor_samples = sorted(donor_samples)

    concordance = {}
    for cs in cluster_samples:
        concordance[cs] = {}
        for ds in donor_samples:
            n_match = 0
            n_compared = 0
            for pos in shared_positions:
                cgt = cluster_gts[pos].get(cs)
                dgt = donor_gts[pos].get(ds)
                if cgt is not None and dgt is not None:
                    n_compared += 1
                    if cgt == dgt:
                        n_match += 1
            concordance[cs][ds] = (n_match, n_compared)

    return concordance, cluster_samples, donor_samples, len(shared_positions)


def greedy_assign(concordance, cluster_samples, donor_samples):
    """Greedy best-first matching: highest concordance pair assigned first."""
    pairs = []
    for cs in cluster_samples:
        for ds in donor_samples:
            n_match, n_compared = concordance[cs][ds]
            rate = n_match / n_compared if n_compared > 0 else 0.0
            pairs.append((rate, n_compared, cs, ds))

    pairs.sort(key=lambda x: (-x[0], -x[1]))

    assigned_clusters = set()
    assigned_donors = set()
    mapping = {}

    for rate, n_compared, cs, ds in pairs:
        if cs in assigned_clusters or ds in assigned_donors:
            continue
        mapping[cs] = (ds, rate, n_compared)
        assigned_clusters.add(cs)
        assigned_donors.add(ds)

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Post-hoc cluster-to-donor assignment via genotype concordance"
    )
    parser.add_argument("--cluster-vcf", required=True,
                        help="Souporcell cluster_genotypes.vcf")
    parser.add_argument("--donor-vcf", required=True,
                        help="WES per-pool donor genotype VCF (decompressed)")
    parser.add_argument("--clusters-tsv", required=True,
                        help="Souporcell clusters.tsv (rewritten in place)")
    parser.add_argument("--report", required=True,
                        help="Output concordance report TSV")
    args = parser.parse_args()

    print(f"Parsing cluster genotypes: {args.cluster_vcf}", file=sys.stderr)
    cluster_samples, cluster_gts = parse_vcf_genotypes(args.cluster_vcf)
    print(f"  {len(cluster_samples)} clusters, {len(cluster_gts)} sites", file=sys.stderr)

    print(f"Parsing donor genotypes: {args.donor_vcf}", file=sys.stderr)
    donor_samples, donor_gts = parse_vcf_genotypes(args.donor_vcf)
    print(f"  {len(donor_samples)} donors, {len(donor_gts)} sites", file=sys.stderr)

    concordance, c_sorted, d_sorted, n_shared = compute_concordance(cluster_gts, donor_gts)
    print(f"  Shared positions: {n_shared}", file=sys.stderr)

    if n_shared == 0:
        print("ERROR: zero shared positions between cluster and donor VCFs", file=sys.stderr)
        sys.exit(1)

    mapping = greedy_assign(concordance, c_sorted, d_sorted)

    # Write full concordance matrix report
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    with open(args.report, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["cluster", "donor", "n_match", "n_compared",
                         "concordance_rate", "assigned"])
        for cs in c_sorted:
            for ds in d_sorted:
                n_match, n_compared = concordance[cs][ds]
                rate = n_match / n_compared if n_compared > 0 else 0.0
                is_assigned = "YES" if cs in mapping and mapping[cs][0] == ds else ""
                writer.writerow([cs, ds, n_match, n_compared, f"{rate:.6f}", is_assigned])

    print("\nCluster-to-donor assignment:", file=sys.stderr)
    for cs in c_sorted:
        if cs in mapping:
            ds, rate, n = mapping[cs]
            print(f"  cluster {cs} -> {ds}  (concordance={rate:.4f}, n={n})", file=sys.stderr)
        else:
            print(f"  cluster {cs} -> UNMATCHED", file=sys.stderr)

    # Validate: warn if any concordance < 0.8
    for cs, (ds, rate, n) in mapping.items():
        if rate < 0.8:
            print(f"WARNING: low concordance for cluster {cs} -> {ds}: "
                  f"{rate:.4f} ({n} sites). Assignment may be unreliable.",
                  file=sys.stderr)

    # Rewrite clusters.tsv in place
    rows = []
    with open(args.clusters_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            a = row["assignment"]
            if "/" in a:
                parts = a.split("/")
                row["assignment"] = "/".join(
                    mapping[p][0] if p in mapping else p for p in parts
                )
            elif a in mapping:
                row["assignment"] = mapping[a][0]
            rows.append(row)

    with open(args.clusters_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Rewrote {args.clusters_tsv} ({len(rows)} barcodes)", file=sys.stderr)
    print(f"Report: {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()
