#!/usr/bin/env python
"""Rename wgsim-simulated reads to unique, informative IDs and build the truth table.

wgsim read names look like:
    @<ref>_<frag_start>_<frag_end>_<errs1>_<errs2>_<pair_idx>/1
where <frag_start>/<frag_end> are 1-based coordinates *local to the FASTA
record wgsim was given* (i.e. relative to the extracted region, not the whole
chromosome). This script:
  1. Converts those to absolute 1-based chromosome coordinates using the
     known 0-based window start (absolute = window_start_0based + local).
  2. Renames each read pair to "<region_id>_<pair_idx>" (unique across
     regions, matching the spirit of spec section 45's suggested read-ID
     convention) and rewrites R1/R2 FASTQ with the new names.
  3. Writes a truth table (spec section 45) with columns:
     read_id, region_id, chromosome, fragment_start, fragment_end,
     read1_start, read2_start (all absolute, 1-based).

Assumes error-free, mutation-free, indel-free simulation (wgsim -e 0 -r 0 -R 0)
so read1_start == fragment_start and read2_start == fragment_end - read_len + 1
exactly (no drift from introduced variants). Read length is read directly from
the FASTQ records, not assumed.
"""
import argparse
import gzip
import re
import sys

# Note: wgsim prints the trailing pair index in hexadecimal (observed sequence
# ...permalink 8, 9, a, b, ... not decimal), so it's matched as a hex string,
# not \d+.
WGSIM_NAME_RE = re.compile(r"^(?P<ref>.+)_(?P<start>\d+)_(?P<end>\d+)_(?P<e1>\d+):(?P<s1>\d+):(?P<i1>\d+)_(?P<e2>\d+):(?P<s2>\d+):(?P<i2>\d+)_(?P<idx>[0-9a-fA-F]+)$")


def open_maybe_gz(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def read_fastq_records(path):
    with open_maybe_gz(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1-in", required=True)
    ap.add_argument("--r2-in", required=True)
    ap.add_argument("--r1-out", required=True)
    ap.add_argument("--r2-out", required=True)
    ap.add_argument("--region-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--window-start-0based", type=int, required=True)
    ap.add_argument("--truth-out", required=True, help="Output truth TSV (appended if it exists)")
    args = ap.parse_args()

    truth_rows = []
    with open_maybe_gz(args.r1_out, "w") as out1, open_maybe_gz(args.r2_out, "w") as out2:
        for (h1, seq1, p1, q1), (h2, seq2, p2, q2) in zip(
            read_fastq_records(args.r1_in), read_fastq_records(args.r2_in)
        ):
            name1 = h1[1:].split()[0]
            assert name1.endswith("/1"), f"expected /1 suffix, got {name1}"
            base_name = name1[:-2]
            m = WGSIM_NAME_RE.match(base_name)
            if not m:
                raise ValueError(f"Could not parse wgsim read name: {base_name!r}")
            local_start = int(m.group("start"))
            local_end = int(m.group("end"))
            pair_idx = m.group("idx")

            abs_frag_start = args.window_start_0based + local_start
            abs_frag_end = args.window_start_0based + local_end
            read1_start = abs_frag_start
            read2_start = abs_frag_end - len(seq2) + 1

            new_id = f"{args.region_id}_{pair_idx}"
            out1.write(f"@{new_id}/1\n{seq1}\n+\n{q1}\n")
            out2.write(f"@{new_id}/2\n{seq2}\n+\n{q2}\n")

            truth_rows.append((new_id, args.region_id, args.chrom, abs_frag_start, abs_frag_end, read1_start, read2_start))

    write_header = True
    try:
        with open(args.truth_out) as f:
            if f.readline():
                write_header = False
    except FileNotFoundError:
        pass

    with open(args.truth_out, "a") as f:
        if write_header:
            f.write("read_id\tregion_id\tchromosome\tfragment_start\tfragment_end\tread1_start\tread2_start\n")
        for row in truth_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wrote {len(truth_rows)} read-pair truth records for {args.region_id} to {args.truth_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
