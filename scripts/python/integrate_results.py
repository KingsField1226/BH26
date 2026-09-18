#!/usr/bin/env python
"""Merge graph_stats.tsv, indexing_stats.tsv, mapping_stats.tsv into combined_results.tsv (spec section 55)."""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-stats", required=True)
    ap.add_argument("--indexing-stats", required=True)
    ap.add_argument("--mapping-stats", required=True)
    ap.add_argument("--regions", required=True, help="metadata/regions.tsv, to restrict to the representative regions and bring in representative_label")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    graph_stats = pd.read_csv(args.graph_stats, sep="\t")
    indexing_stats = pd.read_csv(args.indexing_stats, sep="\t")
    mapping_stats = pd.read_csv(args.mapping_stats, sep="\t")
    regions = pd.read_csv(args.regions, sep="\t")[["representative_label", "region_id"]]

    merged = regions.merge(graph_stats, on="region_id", how="left")
    merged = merged.merge(
        indexing_stats.drop(columns=["window_size"]),
        on="region_id", how="left",
    )
    merged = merged.merge(
        mapping_stats.drop(columns=["window_size"]),
        on="region_id", how="left",
    )
    merged = merged.sort_values("complexity_score").reset_index(drop=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(merged.to_string(index=False))


if __name__ == "__main__":
    main()
