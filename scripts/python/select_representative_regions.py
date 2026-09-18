#!/usr/bin/env python
"""Compute an operational graph-complexity score and select representative regions.

Spec section 36: complexity_score = mean of z-scored
    log10(node_density), log10(edge_node_ratio), log10(sequence_inflation), log10(max_degree + 1)
This is an operational metric for region selection / exploratory comparison
only -- not proposed as a novel biological measure.

Spec section 37: pick 4 representative regions spanning the complexity range
(lowest, lower-middle, upper-middle, highest) and write metadata/regions.tsv.
"""
import argparse
import numpy as np
import pandas as pd


COMPLEXITY_COMPONENTS = ["node_density", "edge_node_ratio", "sequence_inflation", "max_degree"]


def zscore(s):
    mu = s.mean()
    sigma = s.std(ddof=0)
    if sigma == 0 or np.isnan(sigma):
        return pd.Series(0.0, index=s.index)
    return (s - mu) / sigma


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-stats", required=True, help="results/graph_stats.tsv")
    ap.add_argument("--n-representative", type=int, default=4)
    ap.add_argument("--out-scored", required=True, help="Output: graph_stats.tsv augmented with complexity_score")
    ap.add_argument("--out-regions", required=True, help="Output: metadata/regions.tsv (selected representative regions)")
    args = ap.parse_args()

    df = pd.read_csv(args.graph_stats, sep="\t")

    for col in COMPLEXITY_COMPONENTS:
        log_col = f"log10_{col}"
        offset = 1 if col == "max_degree" else 0
        df[log_col] = np.log10(df[col] + offset)

    z_cols = []
    for col in COMPLEXITY_COMPONENTS:
        log_col = f"log10_{col}"
        z_col = f"z_{log_col}"
        df[z_col] = zscore(df[log_col])
        z_cols.append(z_col)

    df["complexity_score"] = df[z_cols].mean(axis=1)
    df = df.sort_values("complexity_score").reset_index(drop=True)
    df.to_csv(args.out_scored, sep="\t", index=False)

    n = args.n_representative
    n_rows = len(df)
    if n_rows < n:
        raise SystemExit(f"Only {n_rows} candidate regions available, need {n} for representative selection")

    # Evenly spaced ranks across the complexity-sorted candidates: lowest, ..., highest.
    rank_positions = [round(i * (n_rows - 1) / (n - 1)) for i in range(n)]
    # De-duplicate while preserving order, in case of collisions with few candidates.
    seen = set()
    chosen_positions = []
    for p in rank_positions:
        if p not in seen:
            chosen_positions.append(p)
            seen.add(p)
    remaining = [i for i in range(n_rows) if i not in seen]
    while len(chosen_positions) < n and remaining:
        # Fill any collision gaps with the next-most-extreme unused rank.
        chosen_positions.append(remaining.pop(0))
    chosen_positions = sorted(chosen_positions)[:n]

    labels = ["A_lowest", "B_lower_middle", "C_upper_middle", "D_highest"]
    if n != 4:
        labels = [f"rank_{i}" for i in range(n)]

    selected = df.iloc[chosen_positions].copy()
    selected.insert(0, "representative_label", labels[: len(selected)])
    selected.to_csv(args.out_regions, sep="\t", index=False)

    print(f"Wrote scored table ({n_rows} regions) to {args.out_scored}")
    print(f"Selected {len(selected)} representative regions to {args.out_regions}")
    print(selected[["representative_label", "region_id", "chrom", "start", "end", "complexity_score"]].to_string(index=False))


if __name__ == "__main__":
    main()
