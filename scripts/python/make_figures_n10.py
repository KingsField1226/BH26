#!/usr/bin/env python
"""Generate Figures 6-8 (n=10 complexity-vs-outcome trend, rank-sampled from
the 326-window corrected pool) from results/combined_results_n10.tsv.
Adapted from make_figures.py's Figure 2-5 style, with region_id labels
(no representative_label concept at n=10) and a Spearman correlation
annotation since a wider sample makes a simple monotonic-trend check useful.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


def scatter_with_labels(ax, df, y_col, y_label, title, highlight=None):
    colors = ["#c0392b" if highlight and r in highlight else "#2b6cb0" for r in df["region_id"]]
    ax.scatter(df["complexity_score"], df[y_col], s=80, color=colors, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(
            row["region_id"].replace("cand_chr21_", "#"),
            (row["complexity_score"], row[y_col]),
            textcoords="offset points", xytext=(6, 4), fontsize=7.5,
        )
    valid = df.dropna(subset=[y_col])
    if len(valid) >= 3:
        rho, pval = spearmanr(valid["complexity_score"], valid[y_col])
        ax.text(0.02, 0.02, f"Spearman rho = {rho:.2f} (p = {pval:.2f})",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.set_xlabel("Complexity score (operational, z-scored; spec section 36)")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined-results", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--highlight", nargs="*", default=[],
                     help="region_ids to highlight in red (known outliers/failures)")
    args = ap.parse_args()

    df = pd.read_csv(args.combined_results, sep="\t").sort_values("complexity_score")

    n_note = (
        f"n = {len(df)} regions rank-sampled evenly across complexity_score from the 326-window "
        "corrected pool (chr21); still descriptive / exploratory, not a statistically powered study."
    )

    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_with_labels(ax, df, "max_rss_gb", "Giraffe index peak RSS (GB)",
                         "Figure 7: Graph complexity vs indexing peak memory (n=10)", args.highlight)
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure7_n10_complexity_vs_index_memory.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_with_labels(ax, df, "wall_time_sec", "Giraffe index wall-clock time (s)",
                         "Figure 6: Graph complexity vs indexing wall time (n=10)", args.highlight)
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure6_n10_complexity_vs_index_time.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_with_labels(ax, df, "accuracy", "Mapping accuracy (correct / mapped)",
                         "Figure 8: Graph complexity vs mapping accuracy (n=10)", args.highlight)
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure8_n10_complexity_vs_accuracy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_with_labels(ax, df, "high_mapq_fraction", "Fraction of mapped reads with MAPQ >= 60",
                         "Figure 9: Graph complexity vs high-MAPQ fraction (n=10)", args.highlight)
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure9_n10_complexity_vs_high_mapq.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote figures 6-9 (n=10) to {args.out_dir}/")


if __name__ == "__main__":
    main()
