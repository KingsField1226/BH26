#!/usr/bin/env python
"""Generate Figure 2: the 3-panel n=50 result figure --
(a) complexity_score vs. mapping accuracy, colored by hazard_flag, with
    pipeline failures shown as markers at the bottom of the plot;
(b) pipeline failure rate by complexity_tier;
(c) mapping accuracy (mean/min) by hazard_flag among successful regions.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def panel_letter(ax, letter):
    ax.text(-0.04, 1.06, letter, transform=ax.transAxes, fontsize=16,
             fontweight="bold", va="top")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined-results", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.combined_results, sep="\t")
    success = df[df["pipeline_status"] == "success"].copy()
    failed = df[df["pipeline_status"] != "success"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    ax = axes[0]
    for hazard, color, label in [
        (False, "#2b6cb0", "No known hazard"),
        (True, "#c0392b", "Known SD / tandem-repeat hazard"),
    ]:
        sub = success[success["hazard_flag_v3"] == hazard]
        ax.scatter(sub["complexity_score"], sub["accuracy"], s=60, color=color,
                   label=label, zorder=3)
    if len(failed):
        ax.scatter(failed["complexity_score"], [0.9800] * len(failed), s=90,
                   color="black", marker="v", zorder=4,
                   label=f"Pipeline failure (n = {len(failed)})")
    ax.set_yticks([0.980, 0.985, 0.990, 0.995, 1.000])
    ax.set_ylim(0.9795, 1.0045)
    ax.set_xlabel("Complexity score (z-scored)")
    ax.set_ylabel("Mapping accuracy")
    ax.set_title("Complexity vs. accuracy (n = 50)", fontsize=12.5)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    panel_letter(ax, "a")

    ax = axes[1]
    fail_rate = df.groupby("complexity_tier").apply(
        lambda g: (g["pipeline_status"] != "success").mean()
    ).reindex(["low", "high"])
    bars = ax.bar(["Low tier", "High tier"], fail_rate.values * 100,
                    color=["#2b6cb0", "#c0392b"])
    for b, v in zip(bars, fail_rate.values * 100):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.5, f"{v:.0f}%",
                 ha="center", fontsize=10)
    ax.set_ylabel("Pipeline failure rate (%)")
    ax.set_title("Failure rate by complexity tier", fontsize=12.5)
    ax.set_ylim(0, max(20, fail_rate.max() * 100 * 1.3))
    ax.grid(True, axis="y", alpha=0.3)
    panel_letter(ax, "b")

    ax = axes[2]
    stats = success.groupby("hazard_flag_v3")["accuracy"].agg(["mean", "min"])
    stats = stats.reindex([False, True])
    x = range(len(stats))
    width = 0.35
    ax.bar([i - width / 2 for i in x], stats["mean"] * 100, width, label="Mean",
            color="#4a7fb5")
    ax.bar([i + width / 2 for i in x], stats["min"] * 100, width, label="Minimum",
            color="#c0392b")
    ax.set_xticks(list(x))
    ax.set_xticklabels(["No known hazard", "Hazard-flagged"])
    ax.set_ylabel("Mapping accuracy (%)")
    ax.set_ylim(97, 100.5)
    ax.set_title("Accuracy by hazard status (successes only)", fontsize=12.5)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    panel_letter(ax, "c")

    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
