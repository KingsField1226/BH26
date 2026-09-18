#!/usr/bin/env python
"""Generate Figures 2-5 (spec section 59) from results/combined_results.tsv.

Figure 1 (low vs high complexity graph topology) and Figure 6 (100kb vs 250kb
comparison, optional) are produced separately (odgi viz / a dedicated script)
since they are not simple scatter plots of this table.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def scatter_with_labels(ax, df, y_col, y_label, title):
    df = df.reset_index(drop=True)
    ax.scatter(df["complexity_score"], df[y_col], s=80, color="#2b6cb0", zorder=3)

    ymin, ymax = df[y_col].min(), df[y_col].max()
    yrange = ymax - ymin
    pad = max(yrange * 0.35, 1e-6)
    ax.set_ylim(ymin - pad, ymax + pad)

    # With only 4 points, the most robust way to avoid label collisions is a
    # vertical zigzag keyed on rank-by-y (not by x/complexity_score order),
    # since points close together in y -- regardless of their x separation --
    # are the ones whose labels actually risk overlapping.
    y_rank = df[y_col].rank(method="first").astype(int) - 1
    dy_pattern = [16, -34, 16, -34]
    for i, row in df.iterrows():
        dy = dy_pattern[y_rank[i] % len(dy_pattern)]
        ax.annotate(
            f"{row['representative_label']}\n({row['region_id']})",
            (row["complexity_score"], row[y_col]),
            textcoords="offset points", xytext=(6, dy), fontsize=8, ha="left",
        )
    ax.set_xlabel("Complexity score (z-scored)")
    ax.set_ylabel(y_label)
    ax.set_title(title, fontsize=12.5)
    ax.grid(True, alpha=0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined-results", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.combined_results, sep="\t").sort_values("complexity_score")

    n_note = (
        f"n = {len(df)} representative 100 kb regions from human chromosome 21; "
        "descriptive comparison, not a statistically powered study."
    )

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter_with_labels(ax, df, "max_rss_gb", "Giraffe index peak memory (RSS, GB)",
                         "Complexity score vs. indexing peak memory")
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure2_complexity_vs_index_memory.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter_with_labels(ax, df, "wall_time_sec", "Giraffe index wall-clock time (s)",
                         "Complexity score vs. indexing wall-clock time")
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure3_complexity_vs_index_time.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter_with_labels(ax, df, "accuracy", "Mapping accuracy (correct / mapped)",
                         "Complexity score vs. mapping accuracy")
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure4_complexity_vs_accuracy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter_with_labels(ax, df, "high_mapq_fraction", "Fraction of mapped reads with MAPQ >= 60",
                         "Complexity score vs. high-confidence (MAPQ >= 60) mapping fraction")
    fig.text(0.5, -0.02, n_note, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/figure5_complexity_vs_high_mapq_fraction.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote figures 2-5 to {args.out_dir}/")


if __name__ == "__main__":
    main()
