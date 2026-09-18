#!/usr/bin/env python
"""Truth-based evaluation of Giraffe mapping (post-surject BAM vs simulated-read truth).

Spec sections 50-51: a mapped read is "correct" if
    |mapped_start - true_start| <= 10 bp
(same chromosome), evaluated per read (both mates of a pair where present).

Input:
  --bam       A BAM/SAM file surjected onto the GRCh38 reference path (samtools view -h compatible).
  --truth     TSV with columns: read_id, region_id, chromosome, fragment_start, fragment_end,
              read1_start, read2_start (see scripts/python/process_wgsim_output.py for how this
              is produced).

Note on correctness check: wgsim does not record which strand a fragment was
sampled from, so whether "read1" (the /1 mate) ends up at the fragment's
start (forward-strand fragment) or at the fragment's end (reverse-strand
fragment) is not knowable from the read name alone -- assuming a fixed
mate/position assignment was tried first and produced ~50% "accuracy" by
construction (exactly the fraction of reverse-strand fragments). Instead,
a mapped mate is scored correct if its absolute position matches *either*
candidate fragment-boundary position (fragment_start, or fragment_end -
read_length + 1) within tolerance, regardless of which mate (read1/read2)
it is. With fragments (~350bp) much longer than reads (150bp) plus the
10bp tolerance, the two candidate positions are never close enough to be
ambiguous in practice.
  --out       Output mapping_stats.tsv row(s) -- appended if the file exists, created otherwise.
  --region-id, --window-size, --complexity-score : identifying columns to attach to the output row
              (complexity_score may be looked up later during results integration instead; pass
              0/NA here if not yet known and patch in during scripts/python/integrate_results.py).

Does not compute wall-clock time or peak memory -- those come from the /usr/bin/time -l wrapper
around the mapping shell script and are merged in separately.
"""
import argparse
import re
import statistics
import subprocess
import sys
import pandas as pd

# odgi extract renames a reference path that only partially covers the
# original chromosome path into a PanSN "subpath" name with an absolute
# 0-based coordinate range suffix, e.g. "GRCh38#0#chr21:38134180-38234292"
# (sometimes with a further "#<fragment>" suffix after that). vg surject's
# BAM RNAME is this full subpath name, and POS is 1-based *relative to the
# start of that subpath*, not an absolute chromosome coordinate. This must
# be converted back to an absolute coordinate before comparing to truth.
SUBPATH_RANGE_RE = re.compile(r"#(?P<contig>[^#:]+):(?P<start>\d+)-(?P<end>\d+)")


def resolve_absolute_position(rname, pos):
    """Return (contig, absolute_1based_pos) from a BAM RNAME + relative POS."""
    m = SUBPATH_RANGE_RE.search(rname)
    if not m:
        # Not a coordinate-range subpath name -- treat rname itself as the
        # contig and pos as already absolute (e.g. a plain "chr21" path).
        return rname, pos
    contig = m.group("contig")
    start_0based = int(m.group("start"))
    return contig, start_0based + pos


def strip_mate_suffix(name):
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def load_truth(path):
    df = pd.read_csv(path, sep="\t", dtype={"read_id": str})
    truth = {}
    for _, row in df.iterrows():
        rid = row["read_id"]
        truth[rid] = {
            "chromosome": str(row["chromosome"]),
            "fragment_start": int(row["fragment_start"]),
            "fragment_end": int(row["fragment_end"]),
        }
    return truth


def iter_sam_records(bam_path):
    """Yield SAM fields for primary alignments only (skip secondary/supplementary)."""
    proc = subprocess.Popen(
        ["samtools", "view", bam_path],
        stdout=subprocess.PIPE, text=True,
    )
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        qname, flag, rname, pos, mapq = fields[0], int(fields[1]), fields[2], int(fields[3]), int(fields[4])
        if flag & 0x100 or flag & 0x800:
            continue  # secondary / supplementary
        yield qname, flag, rname, pos, mapq
    proc.wait()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--region-id", required=True)
    ap.add_argument("--window-size", type=int, required=True)
    ap.add_argument("--tolerance-bp", type=int, default=10)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--high-mapq-threshold", type=int, default=60)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    truth = load_truth(args.truth)

    total_reads = 0
    mapped_reads = 0
    correct_reads = 0
    incorrect_reads = 0
    mapqs = []

    unmatched_truth = 0
    for qname, flag, rname, pos, mapq in iter_sam_records(args.bam):
        base_id = strip_mate_suffix(qname)
        if base_id not in truth:
            unmatched_truth += 1
            continue
        total_reads += 1
        true_chrom = truth[base_id]["chromosome"]
        candidate_positions = (
            truth[base_id]["fragment_start"],
            truth[base_id]["fragment_end"] - args.read_length + 1,
        )

        is_unmapped = bool(flag & 0x4)
        if is_unmapped:
            continue
        mapped_reads += 1
        mapqs.append(mapq)

        contig, absolute_pos = resolve_absolute_position(rname, pos)
        if contig == true_chrom and any(
            abs(absolute_pos - cand) <= args.tolerance_bp for cand in candidate_positions
        ):
            correct_reads += 1
        else:
            incorrect_reads += 1

    if unmatched_truth:
        print(f"WARNING: {unmatched_truth} BAM records had no matching truth entry", file=sys.stderr)

    mapping_rate = mapped_reads / total_reads if total_reads else float("nan")
    accuracy = correct_reads / mapped_reads if mapped_reads else float("nan")
    sensitivity = correct_reads / total_reads if total_reads else float("nan")
    mean_mapq = statistics.mean(mapqs) if mapqs else float("nan")
    median_mapq = statistics.median(mapqs) if mapqs else float("nan")
    high_mapq_fraction = (sum(1 for m in mapqs if m >= args.high_mapq_threshold) / len(mapqs)) if mapqs else float("nan")

    row = {
        "region_id": args.region_id,
        "window_size": args.window_size,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "mapping_rate": mapping_rate,
        "correct_reads": correct_reads,
        "incorrect_reads": incorrect_reads,
        "accuracy": accuracy,
        "sensitivity": sensitivity,
        "mean_mapq": mean_mapq,
        "median_mapq": median_mapq,
        "high_mapq_fraction": high_mapq_fraction,
    }

    out_df = pd.DataFrame([row])
    try:
        existing = pd.read_csv(args.out, sep="\t")
        existing = existing[existing["region_id"] != args.region_id]
        out_df = pd.concat([existing, out_df], ignore_index=True)
    except FileNotFoundError:
        pass
    out_df.to_csv(args.out, sep="\t", index=False)
    print(row, file=sys.stderr)


if __name__ == "__main__":
    main()
