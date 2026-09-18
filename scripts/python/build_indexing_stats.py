#!/usr/bin/env python
"""Parse build_giraffe_index.sh logs into results/indexing_stats.tsv (spec section 41)."""
import argparse
import re
import sys
import pandas as pd


def parse_time_l(log_text):
    """Extract wall-clock seconds and peak RSS (bytes) from `/usr/bin/time -l` output."""
    real_match = re.search(r"^\s*([\d.]+)\s+real", log_text, re.MULTILINE)
    rss_match = re.search(r"^\s*(\d+)\s+maximum resident set size", log_text, re.MULTILINE)
    wall_time_sec = float(real_match.group(1)) if real_match else float("nan")
    max_rss_bytes = int(rss_match.group(1)) if rss_match else None
    return wall_time_sec, max_rss_bytes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regions", required=True, help="metadata/regions.tsv (region_id, window_size, complexity_score, ...)")
    ap.add_argument("--index-dir", required=True, help="work/indexes/")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    regions = pd.read_csv(args.regions, sep="\t")

    rows = []
    for _, region in regions.iterrows():
        region_id = region["region_id"]
        log_path = f"{args.index_dir}/{region_id}/index.log"
        try:
            with open(log_path) as f:
                log_text = f.read()
        except FileNotFoundError:
            rows.append({
                "region_id": region_id, "window_size": region["window_size"],
                "complexity_score": region["complexity_score"],
                "index_success": False, "wall_time_sec": None, "max_rss_gb": None,
                "gbz_size_mb": None, "dist_size_mb": None, "min_size_mb": None,
                "exit_code": None, "failure_reason": "index.log not found",
            })
            continue

        wall_time_sec, max_rss_bytes = parse_time_l(log_text)
        exit_match = re.search(r"INDEX_EXIT_CODE=(-?\d+)", log_text)
        exit_code = int(exit_match.group(1)) if exit_match else None
        index_success = exit_code == 0

        import os
        def size_mb(path):
            try:
                return os.path.getsize(path) / 1e6
            except OSError:
                return None

        gbz_size_mb = size_mb(f"{args.index_dir}/{region_id}/{region_id}.giraffe.gbz")
        dist_size_mb = size_mb(f"{args.index_dir}/{region_id}/{region_id}.dist")
        min_size_mb = size_mb(f"{args.index_dir}/{region_id}/{region_id}.shortread.withzip.min")

        failure_reason = "" if index_success else "see index.log"

        rows.append({
            "region_id": region_id,
            "window_size": region["window_size"],
            "complexity_score": region["complexity_score"],
            "index_success": index_success,
            "wall_time_sec": wall_time_sec,
            "max_rss_gb": (max_rss_bytes / 1e9) if max_rss_bytes else None,
            "gbz_size_mb": gbz_size_mb,
            "dist_size_mb": dist_size_mb,
            "min_size_mb": min_size_mb,
            "exit_code": exit_code,
            "failure_reason": failure_reason,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(out_df.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
