#!/usr/bin/env python
"""Collect per-region odgi stats/degree output into results/graph_stats.tsv.

Reads a candidate-region manifest (region_id, chrom, start, end, window_size,
...) plus, for each region_id, <extract_dir>/<region_id>/odgi_stats.tsv and
odgi_degree.tsv (produced by scripts/shell/extract_region.sh), and writes one
combined table with the raw metrics (spec section 34) plus derived metrics
(spec section 35): node_density, edge_node_ratio, sequence_inflation.
"""
import argparse
import sys
import pandas as pd


def parse_odgi_stats(path):
    """odgi stats -S output: '#length\\tnodes\\tedges\\tpaths\\tsteps' header + one data row."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    row = df.iloc[0]
    return {
        "graph_bp": int(row["length"]),
        "node_count": int(row["nodes"]),
        "edge_count": int(row["edges"]),
        "path_count": int(row["paths"]),
        "step_count": int(row["steps"]),
    }


def parse_odgi_degree(path):
    """odgi degree -S output: '#node.count\\tedge.count\\tavg.degree\\tmin.degree\\tmax.degree'."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    row = df.iloc[0]
    return {
        "mean_degree": float(row["avg.degree"]),
        "max_degree": float(row["max.degree"]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="candidate_regions.tsv (region_id, chrom, start, end, window_size, ...)")
    ap.add_argument("--extract-dir", required=True, help="Directory containing <region_id>/odgi_stats.tsv + odgi_degree.tsv")
    ap.add_argument("--out", required=True, help="Output graph_stats.tsv path")
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")

    rows = []
    for _, region in manifest.iterrows():
        region_id = region["region_id"]
        region_dir = f"{args.extract_dir}/{region_id}"
        stats_path = f"{region_dir}/odgi_stats.tsv"
        degree_path = f"{region_dir}/odgi_degree.tsv"
        try:
            stats = parse_odgi_stats(stats_path)
            degree = parse_odgi_degree(degree_path)
        except FileNotFoundError as e:
            print(f"WARNING: skipping {region_id}, missing output: {e}", file=sys.stderr)
            continue

        reference_length = int(region["end"]) - int(region["start"])
        row = {
            "region_id": region_id,
            "chrom": region["chrom"],
            "start": int(region["start"]),
            "end": int(region["end"]),
            "window_size": int(region["window_size"]),
            "reference_length": reference_length,
            **stats,
            **degree,
        }
        # Derived metrics (spec section 35)
        row["node_density"] = row["node_count"] / row["reference_length"]
        row["edge_node_ratio"] = row["edge_count"] / row["node_count"]
        row["sequence_inflation"] = row["graph_bp"] / row["reference_length"]
        rows.append(row)

    out_df = pd.DataFrame(rows)
    col_order = [
        "region_id", "chrom", "start", "end", "window_size", "reference_length",
        "graph_bp", "node_count", "edge_count", "path_count", "step_count",
        "mean_degree", "max_degree",
        "node_density", "edge_node_ratio", "sequence_inflation",
    ]
    out_df = out_df[[c for c in col_order if c in out_df.columns]]
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(out_df)} regions to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
