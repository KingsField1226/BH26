#!/usr/bin/env python
"""Merge graph_stats.tsv, indexing_stats.tsv, mapping_stats.tsv and the
hazard annotation into the final n=50 combined_results table -- the single
clean-pipeline counterpart of scripts/python/integrate_results.py (which
requires a representative_label column that a stratified sample doesn't have).
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="region_id, chrom, start, end for the 50 regions")
    ap.add_argument("--graph-stats", required=True)
    ap.add_argument("--indexing-stats", required=True)
    ap.add_argument("--mapping-stats", required=True)
    ap.add_argument("--hazard-annotation", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")[["region_id"]]
    graph_stats = pd.read_csv(args.graph_stats, sep="\t")
    indexing_stats = pd.read_csv(args.indexing_stats, sep="\t")
    mapping_stats = pd.read_csv(args.mapping_stats, sep="\t")
    hazard = pd.read_csv(args.hazard_annotation, sep="\t")[[
        "region_id", "n_segdup", "max_segdup_identity",
        "max_tandem_repeat_span_bp", "hazard_flag_v3", "complexity_tier",
    ]]

    merged = manifest.merge(graph_stats, on="region_id", how="left")
    merged = merged.merge(
        indexing_stats.drop(columns=["window_size", "complexity_score"]),
        on="region_id", how="left",
    )
    merged = merged.merge(
        mapping_stats.drop(columns=["window_size"]),
        on="region_id", how="left",
    )
    merged = merged.merge(hazard, on="region_id", how="left")

    merged["pipeline_status"] = merged.apply(
        lambda r: "success" if r.get("index_success") and pd.notna(r.get("accuracy")) else "timeout_failed",
        axis=1,
    )

    merged = merged.sort_values("complexity_score").reset_index(drop=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} rows to {args.out}")
    print(merged[["region_id", "complexity_score", "complexity_tier", "hazard_flag_v3", "pipeline_status", "accuracy"]].to_string(index=False))


if __name__ == "__main__":
    main()
