#!/usr/bin/env python
"""Screen a chromosome-level GRCh38 reference FASTA for candidate windows.

Slides non-overlapping windows of a given size across the chromosome,
computes the N-fraction of each window from the reference FASTA (independent
of the graph), excludes chromosome-end regions and N-rich windows, and
selects a set of candidates spread across the remaining genome.

Output: metadata/candidate_regions.tsv with columns:
    region_id, chrom, start, end, window_size, n_fraction
(path_name is filled in separately once the graph's actual GRCh38 path name
is confirmed from the GFA/ODGI paths listing -- see spec section 28/30.)
"""
import argparse
import sys


def read_fasta_single(path):
    """Read a single-record FASTA file into (header, sequence) with sequence uppercased."""
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    raise ValueError("Expected a single-record FASTA, found a second header")
                header = line[1:]
            else:
                seq_chunks.append(line.upper())
    if header is None:
        raise ValueError("No FASTA header found")
    return header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="Chromosome reference FASTA (single record)")
    ap.add_argument("--chrom", required=True, help="Chromosome label to record in output (e.g. chr21)")
    ap.add_argument("--window-size", type=int, default=100_000)
    ap.add_argument("--end-exclude-fraction", type=float, default=0.02,
                     help="Fraction of chromosome length to exclude from each end")
    ap.add_argument("--max-n-fraction", type=float, default=0.05)
    ap.add_argument("--n-candidates", type=int, default=10,
                     help="Number of evenly-spaced candidates to sample (ignored if --all is set)")
    ap.add_argument("--all", action="store_true",
                     help="Output every window that passes the N-fraction filter, instead of "
                          "sampling --n-candidates evenly-spaced ones. Region IDs are numbered "
                          "by genomic order across all passing windows.")
    ap.add_argument("--out", required=True, help="Output candidate_regions.tsv path")
    args = ap.parse_args()

    header, seq = read_fasta_single(args.fasta)
    chrom_len = len(seq)
    print(f"Loaded {header}: {chrom_len:,} bp", file=sys.stderr)

    exclude_bp = int(chrom_len * args.end_exclude_fraction)
    lo = exclude_bp
    hi = chrom_len - exclude_bp - args.window_size

    windows = []
    for start in range(lo, hi, args.window_size):
        end = start + args.window_size
        window_seq = seq[start:end]
        n_count = window_seq.count("N")
        n_fraction = n_count / len(window_seq)
        windows.append((start, end, n_fraction))

    passing = [w for w in windows if w[2] <= args.max_n_fraction]
    print(f"Total windows: {len(windows)}, passing N<={args.max_n_fraction}: {len(passing)}", file=sys.stderr)

    if args.all:
        chosen_idx = list(range(len(passing)))
    else:
        if len(passing) < args.n_candidates:
            raise SystemExit(f"Only {len(passing)} windows pass the N-fraction filter, "
                              f"need {args.n_candidates}. Loosen --max-n-fraction or --end-exclude-fraction.")

        # Spread candidates evenly across the passing windows (by index order == genomic order).
        step = len(passing) / args.n_candidates
        chosen_idx = sorted(set(int(i * step) for i in range(args.n_candidates)))
        # In case rounding collapsed indices, top up from the remaining pool.
        remaining = [i for i in range(len(passing)) if i not in chosen_idx]
        while len(chosen_idx) < args.n_candidates and remaining:
            chosen_idx.append(remaining.pop(len(remaining) // 2))
        chosen_idx = sorted(chosen_idx)[: args.n_candidates]

    n_digits = max(2, len(str(len(chosen_idx) - 1)))
    with open(args.out, "w") as out:
        out.write("region_id\tchrom\tstart\tend\twindow_size\tn_fraction\n")
        for rank, idx in enumerate(chosen_idx):
            start, end, n_fraction = passing[idx]
            region_id = f"cand_{args.chrom}_{rank:0{n_digits}d}"
            out.write(f"{region_id}\t{args.chrom}\t{start}\t{end}\t{args.window_size}\t{n_fraction:.6f}\n")

    print(f"Wrote {len(chosen_idx)} candidate windows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
