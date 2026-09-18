#!/usr/bin/env python
"""Generate Figure 1: distribution of complexity_score across all 326 valid
chr21 windows, with the actual n=50 stratified sample (the final analysis)
marked. Reads results/graph_stats.tsv (population) and
results/combined_results_n50_stratified.tsv (sample).
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--graph-stats", required=True)
    ap.add_argument("--n50-sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pop = pd.read_csv(args.graph_stats, sep="\t")
    cs = pop["complexity_score"]
    sample = pd.read_csv(args.n50_sample, sep="\t")["complexity_score"]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(cs, bins=30, color="#4a7fb5", edgecolor="white", zorder=2)
    for i, x in enumerate(sample):
        ax.axvline(x, color="#c0392b", linestyle="--", linewidth=0.9,
                    alpha=0.6, zorder=3,
                    label="n = 50 stratified sample (this analysis)" if i == 0 else None)

    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.set_xlabel("Complexity score (z-scored)")
    ax.set_ylabel(f"Number of 100 kb windows (n = {len(cs)})")
    ax.set_title("Distribution of local graph complexity across chr21 windows",
                  fontsize=13)
    ax.grid(True, alpha=0.3, zorder=0)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
