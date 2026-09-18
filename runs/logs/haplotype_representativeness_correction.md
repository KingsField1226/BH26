# Correction: haplotype-representativeness filter, and what it changed

## How this was found

While inspecting `cand_chr21_017.gfa` (the original `A_lowest`) in Bandage, the
user noticed the graph looked sparse in an unexpected way. Counting paths
directly in the GFA showed only **19 of 466** possible haplotype paths were
actually present. Only 2 of the 3 previously-degenerate exclusion criteria
(`edge_count==0`, `path_count<400`, `graph_bp>3x window_size`) had any chance
of catching this, and none did: `cand_chr21_017`'s *screening-stage*
(`-d 0`, unmerged) `path_count` was 1554 -- comfortably above the 400
threshold -- because that count measures fragmented path *pieces*, not
distinct haplotypes. A single haplotype's local assembly can contribute many
small fragments to the unmerged screening graph even when most *other*
haplotypes contribute none at all (a real assembly gap). The filter was
checking the wrong thing.

## Quantifying it properly

For each of the 345 previously-valid windows, distinct haplotypes were
counted directly from the existing screening-stage `.og` files (cheap: just
`odgi paths -L`, parsing the `sample#hap` prefix of each path name, no
re-extraction needed):

```
n_distinct_haplotypes = number of unique "sample#hap" prefixes among all
path names in the window's screening graph
```

Distribution across the 345 windows: median and 75th percentile were both
**466** (the graph's full haplotype count), but 19 windows fell far below
that, with a clean natural gap between 433 and 463:

| range of n_distinct_haplotypes | count of windows |
|---|---|
| 16-22 | 6 |
| 97-128 | 5 |
| 398-433 | 8 |
| 463-466 | 326 |

All 19 low-representation windows are **contiguous on the reference**,
spanning chr21:5,234,199-10,734,199 -- immediately flanking the centromere
(chr21 centromere: 10,864,560-12,915,808 per UCSC). This is a real,
biologically expected pattern: pericentromeric sequence is difficult to
assemble, so a meaningful fraction of the HPRC year-2 draft haplotype
assemblies simply lack a contig spanning parts of this ~5.5 Mb corridor.
`cand_chr21_010` (the rDNA/NOR region, see
`cand_chr21_010_index_failure.md`) sits at the near edge of this same
corridor (432/466 haplotypes) -- a related but separate manifestation of the
same pericentromeric-proximity effect.

## The fix

A new filter, `n_distinct_haplotypes >= 450`, was applied on top of the
existing degenerate-window filter, shrinking the valid candidate pool from
345 to **326** windows (full exclusion list:
`results/graph_stats_excluded_low_haplotype_rep.tsv`; full counts for all
345:  `results/haplotype_representation_all345.tsv`). Complexity scoring and
representative-region selection (`select_representative_regions.py`) were
rerun on this corrected pool.

## What changed

| label | old region_id | old score | new region_id | new score |
|---|---|---|---|---|
| A_lowest | cand_chr21_017 | -1.099 | **cand_chr21_050** | **-2.166** |
| B_lower_middle | cand_chr21_346 | -0.251 | **cand_chr21_093** | **-0.236** |
| C_upper_middle | cand_chr21_274 | 0.044 | **cand_chr21_287** | **0.109** |
| D_highest | cand_chr21_330 | 0.630 | cand_chr21_330 (unchanged) | 0.630 |

**D_highest was deliberately NOT changed.** The pool's new strict top pick is
`cand_chr21_347` (score 2.889) -- but this is already a documented failure
from the `D_highest` substitution saga (killed at the `vg giraffe` mapping
stage; see `D_highest_substitution_saga.md`). Rather than re-launch another
multi-attempt substitution search, `cand_chr21_330` (proven working twice
over, and still at the 90.8th percentile -- rank 30/326 -- of the corrected
pool) was kept. For the record, the new pool's top 10 by score includes at
least 3 previously-confirmed failures (`cand_chr21_347` rank 1,
`cand_chr21_354` rank 5, `cand_chr21_117` rank 10) -- further evidence that
high complexity_score and pipeline failure are correlated but that chasing
the single highest score is not a reliable strategy at this scale.

## A second problem found along the way: `cand_chr21_093` also blew up

`cand_chr21_093`'s *final* extraction (default `-d 300000` subpath merging)
ballooned to **537,760 bp -- 5.4x** the requested 100 kb window (screening-
stage size was a normal 100,032 bp), the same "merge chases into an adjacent
collapsed-repeat structure" pathology previously seen for `cand_chr21_010`
(46x) and `cand_chr21_354` (2.5x). Despite the blowup, indexing and mapping
both completed successfully (100% accuracy) -- but a 538 kb "B_lower_middle"
breaks the project's core "same window size, varying complexity" comparison
design. Re-extracted with `merge_distance=0` (same override used for
`cand_chr21_010`) to get the correct 100,032 bp size; the original blown-up
extraction is kept for the record at
`work/extracted_final/cand_chr21_093_blownup_default_merge/`. Indexing and
mapping were rerun on the correctly-sized version (still 100% accuracy;
`path_count` is now 66,946 fragmented pieces rather than ~466 haplotype-
length paths, the known cost of disabling subpath merging).

## A striking new finding: the corrected `A_lowest` sits on a segmental-duplication hotspot

`cand_chr21_050` (chr21:13,834,199-13,934,199) produced the **worst mapping
accuracy of all 4 final regions** (98.695%, vs. >=99.87% for B/C and 99.56%
for D_highest) -- worse even than `D_highest`, despite scoring as the single
*lowest*-complexity region in the corrected pool (-2.166, the most extreme
score of any of the 4). This breaks the otherwise-clean "complexity score
predicts accuracy" pattern seen in every prior iteration of this analysis
(pilot and both full-scan drafts).

Checking annotation explains why: UCSC `genomicSuperDups` returns **41
segmental-duplication entries** overlapping this window, at 90-97% sequence
identity, matching paralogous sequence on **chr13, chr9, chr2, chr20, chr14,
chr4, and chr22** (plus unplaced contigs). This is a genuine, classic
inter-chromosomal segmental-duplication hotspot -- a sequence class well
known to cause short-read mapping ambiguity, because reads from here have
near-identical sequence at multiple genomic locations. None of the 5 prior
`D_highest`-search failures had any segmental duplication at all (checked
directly, see `D_highest_substitution_saga.md`), which makes this the first
segmental-duplication-driven finding in the whole project, and it appears
at the *low*-complexity end rather than the high end.

**Why did the complexity score miss this?** The composite score is an
unweighted mean of 4 z-scored metrics computed from the *local* extracted
subgraph's own topology (`node_density`, `edge_node_ratio`,
`sequence_inflation`, `max_degree`). A segmental duplication's defining
property -- near-identical sequence *elsewhere in the genome* -- is
invisible to a metric computed purely from this one local subgraph; if the
466 chr21 haplotypes themselves happen to track GRCh38 closely through this
particular window (they do: `node_density` here is actually the *highest*
of the 4 final regions at 0.112, but `sequence_inflation` is unusually low
at 0.754, pulling the unweighted mean score down overall), the score reads
"simple" even though the locus is one of the more mapping-hazardous regions
in the genome by a different, more established measure. This is a distinct
failure mode from the `D_highest` saga (which was about *pipeline
tractability*, not accuracy) but the same underlying lesson: this project's
composite complexity score, built from local graph topology alone, does not
capture every axis of genomic difficulty -- segmental duplication being a
clear example missed at the low-complexity end, just as several tandem-
repeat/rDNA/NUMT cases were initially missed or required careful
localization at the high-complexity end.
