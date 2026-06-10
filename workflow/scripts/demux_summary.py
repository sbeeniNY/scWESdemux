"""
Aggregate per-pool Souporcell `clusters.tsv` outputs into one cohort-wide summary.

Invoked by Snakemake (`rules/demux.smk :: rule demux_summary`) which provides:
    snakemake.input.clusters    : list of paths to clusters.tsv (one per pool)
    snakemake.input.config_yaml : path to config/config.yaml
    snakemake.output.tsv        : output TSV path
    snakemake.log[0]            : log path

Souporcell clusters.tsv columns:
    barcode, status, assignment, log_prob_singleton, log_prob_doublet

When --known_genotypes_sample_names is used, the `assignment` column
contains the donor names (e.g. DonorA) rather than cluster indices.

Output columns (compatible with v1/Vireo format):
    pool, donor_id, mapped_label, n_cells, doublet_rate, unassigned_rate, total_cells
"""

import os
import re
import sys
import logging
import yaml
import pandas as pd


log_path = snakemake.log[0] if snakemake.log else None  # noqa: F821
if log_path:
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
else:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

log = logging.getLogger("demux_summary")

config_path = snakemake.input.config_yaml  # noqa: F821
cluster_paths = list(snakemake.input.clusters)  # noqa: F821
out_tsv = snakemake.output.tsv  # noqa: F821

with open(config_path, "r") as fh:
    cfg = yaml.safe_load(fh)

sample_labels = cfg.get("sample_labels", {}) or {}
pool_donors = {p: info["donors"] for p, info in cfg.get("pools", {}).items()}


def _pool_from_path(path: str) -> str:
    parts = os.path.normpath(path).split(os.sep)
    try:
        idx = parts.index("souporcell")
        return parts[idx - 1]
    except ValueError:
        m = re.search(r"/demux/([^/]+)/souporcell/clusters\.tsv$", path)
        if not m:
            raise ValueError(f"Cannot infer pool name from path: {path}")
        return m.group(1)


rows: list[dict] = []

for tsv_path in cluster_paths:
    pool = _pool_from_path(tsv_path)
    log.info(f"Reading {pool}: {tsv_path}")

    df = pd.read_csv(tsv_path, sep="\t")

    if "status" not in df.columns or "assignment" not in df.columns:
        log.error(f"{tsv_path} missing expected columns; columns={list(df.columns)}")
        continue

    total = len(df)

    n_doublet = int((df["status"] == "doublet").sum())
    n_unassigned = int((df["status"] == "unassigned").sum())
    doublet_rate = (n_doublet / total) if total else 0.0
    unassigned_rate = (n_unassigned / total) if total else 0.0

    singlets = df[df["status"] == "singlet"]

    donors = pool_donors.get(pool, [])
    assignment_vals = singlets["assignment"].value_counts(dropna=False).to_dict()

    donor_map = {}
    for val in assignment_vals:
        if pd.isna(val):
            continue
        s = str(val)
        if s.isdigit():
            idx = int(s)
            donor_map[val] = donors[idx] if idx < len(donors) else f"cluster_{idx}"
        else:
            donor_map[val] = s

    pool_label_map = sample_labels.get(pool, {}) or {}

    # --- per-pool clusters_with_samples.tsv ---
    def _map_row(row):
        if row["status"] == "doublet":
            return pd.Series({"sample_assignment": "doublet", "sample_label": ""})
        if row["status"] == "unassigned":
            return pd.Series({"sample_assignment": "unassigned", "sample_label": ""})
        a = str(row["assignment"])
        if a.isdigit():
            idx = int(a)
            pid = donors[idx] if idx < len(donors) else f"cluster_{idx}"
        else:
            pid = a
        return pd.Series({
            "sample_assignment": pid,
            "sample_label": pool_label_map.get(pid, ""),
        })

    annotated = pd.concat([df, df.apply(_map_row, axis=1)], axis=1)
    annotated_path = os.path.join(os.path.dirname(tsv_path), "clusters_with_samples.tsv")
    annotated.to_csv(annotated_path, sep="\t", index=False)
    log.info(f"Wrote annotated clusters for {pool}: {annotated_path}")

    for assignment_val, n_cells in assignment_vals.items():
        if pd.isna(assignment_val):
            continue
        donor_id = donor_map.get(assignment_val, str(assignment_val))
        rows.append(
            {
                "pool": pool,
                "donor_id": donor_id,
                "mapped_label": pool_label_map.get(donor_id, ""),
                "n_cells": int(n_cells),
                "doublet_rate": round(doublet_rate, 6),
                "unassigned_rate": round(unassigned_rate, 6),
                "total_cells": total,
            }
        )

    for special, count in [("doublet", n_doublet), ("unassigned", n_unassigned)]:
        if count > 0:
            rows.append(
                {
                    "pool": pool,
                    "donor_id": special,
                    "mapped_label": "",
                    "n_cells": count,
                    "doublet_rate": round(doublet_rate, 6),
                    "unassigned_rate": round(unassigned_rate, 6),
                    "total_cells": total,
                }
            )

summary = pd.DataFrame(
    rows,
    columns=[
        "pool",
        "donor_id",
        "mapped_label",
        "n_cells",
        "doublet_rate",
        "unassigned_rate",
        "total_cells",
    ],
).sort_values(["pool", "donor_id"]).reset_index(drop=True)

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
summary.to_csv(out_tsv, sep="\t", index=False)

log.info(f"Wrote {len(summary)} rows to {out_tsv}")
print(f"[demux_summary] {len(summary)} rows -> {out_tsv}", file=sys.stderr)
