#!/usr/bin/env python3
"""
Standalone cohort-wide demux summary for Souporcell (no Snakemake).

Reads per-pool Souporcell clusters.tsv files and config.yaml sample_labels,
writes a single TSV matching the Snakemake workflow output.

Example:
  python3 scripts/demux_summary_standalone.py \\
    --config config/config.yaml \\
    --pools Pool1 Pool3 \\
    --input-dir /path/to/demux \\
    --output /path/to/demux_summary.tsv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys

import pandas as pd
import yaml

log = logging.getLogger("demux_summary_standalone")


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


def aggregate(config_path: str, cluster_tsv_paths: list[str], out_tsv: str) -> None:
    with open(config_path, "r") as fh:
        cfg = yaml.safe_load(fh)

    sample_labels = cfg.get("sample_labels", {}) or {}
    pool_donors = {p: info["donors"] for p, info in cfg.get("pools", {}).items()}
    rows: list[dict] = []

    for tsv_path in cluster_tsv_paths:
        pool = _pool_from_path(tsv_path)
        log.info("Reading %s: %s", pool, tsv_path)

        if not os.path.isfile(tsv_path):
            log.warning("Missing clusters.tsv (skip): %s", tsv_path)
            continue

        df = pd.read_csv(tsv_path, sep="\t")

        if "status" not in df.columns or "assignment" not in df.columns:
            log.error("%s missing expected columns; columns=%s", tsv_path, list(df.columns))
            continue

        total = len(df)
        n_doublet = int((df["status"] == "doublet").sum())
        n_unassigned = int((df["status"] == "unassigned").sum())
        doublet_rate = (n_doublet / total) if total else 0.0
        unassigned_rate = (n_unassigned / total) if total else 0.0

        singlets = df[df["status"] == "singlet"]
        assignment_vals = singlets["assignment"].value_counts(dropna=False).to_dict()

        donors = pool_donors.get(pool, [])
        donor_map = {}
        for val in assignment_vals:
            if isinstance(val, (int, float)) and not pd.isna(val):
                idx = int(val)
                donor_map[val] = donors[idx] if idx < len(donors) else f"cluster_{idx}"
            else:
                donor_map[val] = str(val)

        pool_label_map = sample_labels.get(pool, {}) or {}

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

    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    summary.to_csv(out_tsv, sep="\t", index=False)
    log.info("Wrote %s rows to %s", len(summary), out_tsv)
    print(f"[demux_summary_standalone] {len(summary)} rows -> {out_tsv}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build demux_summary.tsv from Souporcell clusters.tsv files.")
    parser.add_argument("--config", required=True, help="Path to config.yaml (for sample_labels + pool donors).")
    parser.add_argument("--pools", nargs="+", required=True, help="Pool names, e.g. Pool1 Pool3 ...")
    parser.add_argument("--input-dir", required=True, help="Demux root containing <pool>/souporcell/clusters.tsv")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    paths = [os.path.join(args.input_dir, p, "souporcell", "clusters.tsv") for p in args.pools]
    aggregate(args.config, paths, args.output)


if __name__ == "__main__":
    main()
