# Exploring the complexity limits of short-read mapping to human PGGB graphs

BioHackathon proof-of-concept analysis. **All computation for this project runs
locally on a single Mac** (Apple Silicon, 16 cores, 128 GB unified memory) —
no cloud, remote, or external compute was used at any point (see
"Technical limitations" and the analysis spec for the rationale).

## Project objective

Investigate how local topological complexity of a human PGGB pangenome graph
affects the computational cost (Giraffe index construction time/memory) and
accuracy (short-read mapping) of short-read mapping with `vg giraffe`.

## Project structure

All analysis lives under `runs/`. The design, run as a
single clean pass: *every* N-filter-passing 100 kb window on chr21 was
screened (367 windows, not a sample); each valid window was scored for local
graph complexity and checked against a known-hazard sequence annotation
(segmental duplications, large tandem repeats); a stratified sample of 50
regions was drawn across the resulting 2x2 design (complexity tier x hazard
status); and the full indexing/mapping/evaluation pipeline was run on all 50.

**Headline finding:** local-topology complexity score reliably predicts
whether the pipeline completes at all (0% failure among the 25
below-median-complexity regions, 16% among the 25 above-median regions),
while the known-hazard annotation reliably predicts mapping-accuracy
degradation among regions that do complete (minimum accuracy 99.53% without
a hazard flag vs. 98.69% with one) — and neither predicts the other's
outcome. See "Results" below.

This design was arrived at after an earlier, smaller round of work (a
4-region hand-picked comparison) turned up two serious problems with
picking regions by hand — see "How this design was arrived at" near the end
of Results for the condensed history, and
`runs/logs/` for the full blow-by-blow narratives this
project kept rather than discarding.

Top-level `data/`, `scripts/`, `software/`, and environment-record assets
are shared across the project; run-specific inputs/outputs/results live
under `runs/` (`metadata/`, `results/`, `logs/` — the
`work/` directories of large intermediate pipeline artifacts are removed
once each stage's final results are captured in `results/*.tsv`).

## Background and significance

Every component used here — PGGB, the vg toolkit, Giraffe, the HPRC
pangenome — is established, published work. What appears to still be open,
based on a (non-exhaustive, September 2026) literature check, is the
specific question this project asks: *does a region's local graph topology
predict Giraffe's indexing cost and mapping accuracy there*, measured with
an explicit, quantitative per-window complexity score across many small
regions of a real HPRC PGGB graph.

A few data points from that check, in order of how directly they bear on
this project:

- **The vg team's own most recent Giraffe paper (Sept 2025) states they have
  not been able to build genome-scale Giraffe indexes for PGGB graphs at
  all**, and benchmark Giraffe exclusively on Minigraph-Cactus graphs
  instead: *"PGGB's reference-free approach to graph building provides an
  unbiased view of the homology within and among genomes. However, due to
  the complex topology of the PGGB graphs, we have not been able to
  construct indexes for Giraffe."* ([Liao et al., 2025, bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.29.678807v1);
  see also the [PMC version](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621775/)).
  This project's central finding — that Giraffe indexing and mapping *are*
  tractable for **local** (100 kb) PGGB extracts most of the time, with a
  measurable local complexity score predicting the fraction that will fail
  to complete at all — sits directly against this backdrop: it is a
  small-scale, local-window probe of exactly the failure mode the tool's
  own developers describe encountering at genome scale, rather than a
  replication of already-published results.
- The **original Giraffe paper** already names "graphs with complex
  topology" as a design target Giraffe's haplotype-constrained search is
  meant to handle efficiently ([Sirén et al., 2021, *Science*](https://pubmed.ncbi.nlm.nih.gov/34914532/)),
  but evaluates this at the level of whole-genome graph design choices, not
  via a per-region complexity score correlated with per-region outcomes.
- The **vg wiki's Index Construction page** documents, as long-standing
  operational knowledge, that "if the graph is more complex than expected,
  \[index\] construction may use orders of magnitude more space before
  ultimately failing" ([vgteam/vg wiki](https://github.com/vgteam/vg/wiki/Index-Construction)) —
  corroborating that graph complexity is a recognized general risk factor
  for indexing, though this documentation concerns the older GCSA index
  rather than Giraffe's GBZ/minimizer/distance index, and offers no
  quantitative, per-region complexity measure.
- Work on **defining "easy" vs. "difficult" regions from pangenome data**
  exists, but by a different route: [Li, 2025 (arXiv 2507.03718)](https://arxiv.org/pdf/2507.03718)
  classifies regions via k-mer uniqueness (a ropebwt3 index) for short-read
  *variant-calling* accuracy, not via graph-topology metrics, and does not
  examine Giraffe indexing cost or MAPQ.
- Work **comparing graph topology across construction methods** (PGGB vs.
  Minigraph-Cactus vs. others) exists at whole-graph and specific-locus
  granularity ([Andreace et al., 2023, *Genome Biology*](https://pubmed.ncbi.nlm.nih.gov/38037131/)),
  and a related sliding-window *graph edit distance between methods*
  reveals local "edit hotspots" ([PMC12187062](https://pmc.ncbi.nlm.nih.gov/articles/PMC12187062/))
  — a similar sliding-window spirit to this project's screening step, but
  measuring divergence *between* construction methods rather than one
  graph's own topology and its relationship to Giraffe's indexing/mapping
  behavior.

None of this rules out narrower prior work using this exact combination of
metrics that a deeper review might surface. But the balance of evidence —
including a direct admission from the tool's own current developers that
genome-scale PGGB+Giraffe indexing is presently unsolved — suggests this
project's specific angle (small, local, quantitative, complexity-vs-outcome)
is a genuine, if modest, contribution rather than a restatement of settled
results.

## Research questions

Primary: *At what level of local graph complexity does short-read mapping to
a human PGGB graph become computationally difficult or inaccurate?*

Deeper question the analysis is designed to speak to: *Is short-read mapping
feasibility determined primarily by graph size, graph topology, or both?*
This project addresses the topology side directly, by comparing regions of
differing local graph complexity at a fixed 100 kb window size.

See the full hypothesis set (H1/H2/H3) and question list (Q1-Q11) in the
project's analysis specification document.

## Computer environment

```
OS: macOS 26.6.2 (build 25G83)
Architecture: arm64 (Apple Silicon)
CPU: 16 cores (standard job thread count: 12)
Memory: 128 GB unified memory (target ceiling for analysis jobs: ~90 GB;
         hard stop if a job's RSS approaches ~100-105 GB and keeps climbing)
```

Verified at analysis start with `scripts/shell/check_environment.sh`
(`uname -m`, `sw_vers`, `$CONDA_DEFAULT_ENV`, `which python`).

## Conda environment

All analysis tools live in a single dedicated conda environment, `bh26`
(Python 3.11, conda-forge + bioconda channels). Nothing was installed via
Homebrew, MacPorts, or system Python; no other conda environment, and no
system/Homebrew configuration, was modified. Full package list and exported
environment specs are under `results/`:
`conda_list_bh26.txt`, `bh26_environment.yml`, `bh26_environment_from_history.yml`.

Key tools: `odgi`, `samtools`, `seqkit`, `wgsim`, `vg` (built from source —
see below), plus Python `pandas`/`numpy`/`matplotlib`/`scipy`/`psutil`.

## Input data

HPRC (Human Pangenome Reference Consortium) **release2** PGGB pangenome
graph, **chromosome-level GFA for chr21** (parameters `p98-k311`), downloaded
from the public `human-pangenomics` S3 bucket. This graph includes GRCh38 and
CHM13 as reference paths plus 464 HPRC year-2 haplotypes (466 total
haplotype/reference paths; see `hpp_pangenome_resources` for the release
description). chr21 was chosen as the smallest autosome (~2.4 GB compressed
GFA) to keep the "chromosome-level" starting point (spec priority 2, since a
usable partition-level-with-known-coordinates alternative was not readily
identifiable — see Technical limitations) tractable for local extraction on
this Mac. This `chr21.og` build underlies all screening and analysis in this
project.

An independent GRCh38 chr21 reference FASTA (UCSC `hg38`) was also downloaded,
used only for N-content screening of candidate windows and as a coordinate
cross-check independent of the graph. Two UCSC annotation tracks
(`genomicSuperDups`, `simpleRepeat`) were also downloaded, once each for the
whole of chr21, for the known-hazard annotation described under Methods.

Full provenance (URL, download date, file size, SHA-256) for every downloaded
file is recorded in `metadata/input_files.tsv`.

## Data provenance

See `metadata/input_files.tsv` (machine-readable) for exact source URLs,
checksums, and download dates of every raw input file.

## Software versions

See `results/software_versions.txt`, generated at analysis start and updated
once `vg` finished building.

## Methods

Pipeline, executed in order:

```
HPRC PGGB chromosome graph (chr21, release2, p98-k311)
    -> ODGI graph (odgi build)
    -> screen every 100 kb window (odgi extract/stats/degree + Python)
    -> local graph complexity score, per window
    -> known-hazard annotation, per window (UCSC genomicSuperDups/simpleRepeat)
    -> stratified sample of 50 regions (complexity tier x hazard status)
    -> local graph extraction (odgi extract, GRCh38#0#chr21 coordinates)
    -> Giraffe index construction (vg autoindex)
    -> simulated paired-end reads (wgsim, error-free, seed=42)
    -> Giraffe mapping (vg giraffe)
    -> surjection to GRCh38 (vg surject)
    -> truth-based mapping evaluation
```

### Candidate window screening

*Every* 100 kb window on chr21 was screened
(`scripts/python/screen_candidate_windows.py --all`), excluding the
first/last 2% of the chromosome and any window with >5% N-content in the
GRCh38 reference FASTA. This produced 448 total 100 kb windows, of which 367
pass the N<=5% filter (`runs/metadata/candidate_regions_all.tsv`).
All 367 were screened with `odgi extract`/`odgi stats`/`odgi degree`
(`scripts/shell/run_screen_batch_safe.sh`, ~6.5 hours wall time), which
wraps each window in a hard per-window wall-clock kill so one pathological
window cannot silently stall the whole scan (see "Technical limitations").

Two filters were applied before scoring, leaving **326 valid windows**
(`runs/results/graph_stats.tsv`) used for sampling:

1. **Degenerate-window filter**: `edge_count==0 OR path_count<400 OR
   graph_bp>3x window_size` excludes 22 windows clustered at chr21
   ~6.0-6.1 Mb, ~7.9-8.6 Mb, and ~10.9-12.8 Mb (the pericentromeric/
   acrocentric region), leaving 345 windows.
2. **Haplotype-representativeness filter**: `path_count<400` above is a
   *screening-stage* (unmerged) path-fragment count, not a distinct-
   haplotype count — a handful of haplotypes with heavily fragmented local
   assemblies can produce many fragments while most other haplotypes are
   simply absent from that window. Counting **distinct haplotypes directly**
   (`odgi paths -L`, parsed by `sample#hap` prefix) across the 345 windows
   found 19 with severely reduced representation (as low as 16/466 present),
   forming a contiguous corridor at chr21:5,234,199-10,734,199 immediately
   flanking the centromere — a real, expected consequence of pericentromeric
   sequence being difficult for draft haplotype assemblies to resolve. A
   filter of `n_distinct_haplotypes >= 450` (a clean natural gap in the data
   between 433 and 463) removes these 19, leaving 326 windows. Full counts:
   `results/haplotype_representation_all345.tsv` and
   `results/graph_stats_excluded_low_haplotype_rep.tsv`.

`odgi extract`'s `-E/--full-range` flag was found to be impractical on this
graph (see "Technical limitations") and is not used anywhere in this
pipeline; screening extraction also disables default subpath-merging
(`-d 0`) for speed, while final region extraction (below) uses the default
merging distance.

### Graph complexity metrics

For each candidate window: `reference_length`, `graph_bp`, `node_count`,
`edge_count`, `path_count`, `step_count`, `mean_degree`, `max_degree` (via
`odgi stats -S` / `odgi degree -S`), plus derived metrics `node_density`,
`edge_node_ratio`, `sequence_inflation`. See `results/graph_stats.tsv`.

An operational (not biological) complexity score was computed as the mean of
z-scored `log10(node_density)`, `log10(edge_node_ratio)`,
`log10(sequence_inflation)`, `log10(max_degree + 1)`, via
`scripts/python/select_representative_regions.py`. (Its `zscore()` helper
has a latent bug — dividing by a `sigma` that becomes `NaN` if any input
column contains `-inf` — worked around by pre-filtering degenerate rows
before scoring rather than patched at the source; see "Technical
limitations".)

### Known-hazard annotation

All 326 valid windows were checked against two UCSC tracks
(`genomicSuperDups`, `simpleRepeat`), downloaded once for the whole of chr21
and intersected locally in pandas rather than queried per window. A binary
`hazard_flag` was defined as *any overlapping segmental duplication at
>=90% identity, or a tandem-repeat array >=1.4 kb* — thresholds calibrated
to correctly flag every problem window this project independently
identified during earlier work (`cand_chr21_050`, `_354`, `_329`, `_358`;
see `runs/logs/hazard_stratified_robustness_analysis.md`
for the calibration detail). This flags 77/326 windows (23.6%);
full annotation in `results/hazard_annotation_326windows.tsv`. Crossed with
`complexity_tier` (above/below the population median `complexity_score`,
median = -0.090), this gives a clean 2x2 design.

### Stratified region selection

50 regions were rank-sampled *within each of the 4*
`{complexity_tier} x {hazard_flag}` *cells* (not across the whole pool), so
that hazard-positive-but-low-complexity and hazard-negative-but-high-
complexity windows — exactly the cases a plain complexity-rank sample would
rarely include — are guaranteed representation (25 regions per tier, 25 per
hazard status by construction). Final picks:
`runs/metadata/manifest_n50_final.tsv`. All 50 regions
were run through the full pipeline in a single clean pass — see Results.

### Indexing procedure

`vg autoindex --workflow sr-giraffe` on each region's GFA. Resource usage
recorded via `/usr/bin/time -l`. A **30-minute per-region wall-clock budget**
is enforced end to end (extraction through evaluation) via
`scripts/shell/run_full_pipeline_final_n50.sh`, a single reusable script
chaining extraction -> indexing -> simulation -> mapping -> surjection ->
evaluation with explicit stage-completion markers in each region's log, so a
failure's exact stage is always identifiable. The script also detects and
automatically retries (with `merge_distance=0`) any extraction whose result
exceeds 3x the requested window size — a recurring pathology on this graph
(see "Technical limitations"). This 30-minute budget, and the automatic
blowup retry, are both design choices carried forward from lessons learned
during earlier development of this methodology (see "Technical
limitations" and the resource-budget caveat under "Results"). See
`results/indexing_stats_n50_final.tsv`.

### Read simulation

`wgsim`, paired-end, 150 bp reads, ~350 bp fragments (`-d 350 -s 35`),
error-free/mutation-free/indel-free (`-e 0 -r 0 -R 0`), fixed seed (42),
30,000 read pairs per region, from the GRCh38 reference interval of each
region. Reads renamed to `<region_id>_<pair_idx>` and a machine-readable
truth table written alongside (`scripts/python/process_wgsim_output.py`).
Note: wgsim's trailing read-pair index in the FASTQ read name is printed in
**hexadecimal** (observed sequence `...8, 9, a, b...`), not decimal as might
be assumed — handled by matching `[0-9a-fA-F]+` rather than `\d+`.

### Mapping procedure

`vg giraffe` against each region's GBZ index, followed by `vg surject` back
onto the GRCh38 reference path so mapped coordinates are directly comparable
to simulation truth coordinates. Because each region's GRCh38 path is itself
a *subpath* of the whole chr21 path (PanSN "subpath" naming, e.g.
`GRCh38#0#chr21:38134180-38234292`), `vg paths -x <gbz> -L` resolves the
exact subpath name at run time — it is not the same string for every region
— and downstream BAM positions are relative to that subpath's own start,
not the chromosome, until converted back (see next section).

### Mapping evaluation

A mapped read is scored "correct" if, after converting its BAM position back
to an absolute GRCh38 coordinate (subpath start offset parsed out of the
RNAME, plus the relative POS), it lands on the correct chromosome within
10 bp of a true fragment-boundary position. Metrics: mapping rate, accuracy,
sensitivity, mean/median MAPQ, high-MAPQ fraction (MAPQ>=60). See
`scripts/python/evaluate_mapping.py` and `results/mapping_stats_n50_final.tsv`.

Correctness note: wgsim does not record which strand a simulated fragment
came from, so whether the "/1" mate ends up at the fragment's start or its
end is not recoverable from the read name alone. A first pass that assumed
a fixed assignment (read1=start, read2=end) produced ~50% "accuracy" by
construction — exactly the fraction of reverse-strand fragments — which was
the first sign something was wrong with the check itself, not the mapper.
The fix: a mapped mate is scored correct if its absolute position matches
*either* candidate fragment-boundary position, regardless of which mate it
is; with a 350 bp fragment and 150 bp reads, the two candidate positions are
always far enough apart (~200 bp) that this introduces no real ambiguity
against the 10 bp tolerance.

## Selected regions

`runs/metadata/manifest_n50_final.tsv` (50 regions,
stratified across complexity tier and hazard status).

## Results

Full merged table: `runs/results/combined_results_n50_final.tsv`
(50 rows: `graph_stats.tsv` + `indexing_stats_n50_final.tsv` +
`mapping_stats_n50_final.tsv` + `hazard_annotation_326windows.tsv`, built by
`scripts/python/build_combined_results_n50_final.py`).

**Pipeline completion, by complexity tier (25 regions per tier):**

| | Low tier (below median) | High tier (above median) |
|---|---|---|
| Failed to complete | 0 / 25 (0%) | 4 / 25 (16%) |

**Pipeline completion, by known-hazard status (25 regions per group):**

| | No known hazard | Hazard-flagged |
|---|---|---|
| Failed to complete | 1 / 25 (4%) | 3 / 25 (12%) |

**Mapping accuracy among the 46 regions that completed, by known-hazard
status:**

| | No known hazard (n=24) | Hazard-flagged (n=22) |
|---|---|---|
| Mean accuracy | 99.92% | 99.84% |
| Minimum accuracy | **99.53%** | **98.69%** |

The minimum-accuracy region overall, `cand_chr21_050`
(complexity_score -2.166, the single lowest-scoring region in the entire
326-window pool), carries a hazard flag: UCSC `genomicSuperDups` shows 41
entries at 90-97% identity matching paralogous sequence on chr13, chr9,
chr2, chr20, chr14, chr4, and chr22 — a genuine inter-chromosomal segmental
duplication invisible to a complexity score computed purely from one local
subgraph. **This is the headline pattern**: `complexity_tier` cleanly
separates *whether the pipeline finishes at all*; `hazard_flag` cleanly
separates *how accurate the result is, given that it finishes*; neither
predicts the other outcome. All 4 pipeline failures fell in the
high-complexity tier regardless of hazard status (one of the 4,
`cand_chr21_313`, carries no hazard flag at all), while the single worst
accuracy result belongs to the *lowest*-complexity region in the dataset.

See Figure 1 (`results/figures/figure1_complexity_score_distribution.png`,
the right-skewed distribution of `complexity_score` across all 326 windows
— 77.6% fall within -0.5 to 0.5 — that motivated a stratified rather than
uniform sample) and Figure 2
(`results/figures/figure2_hazard_vs_complexity_2x2.png`, three panels: the
full n=50 sample colored by hazard status with failures marked, failure
rate by complexity tier, and accuracy by hazard status).

**Two important caveats on the numbers above**, both load-bearing for how
strongly to read them:

- **"Failed to complete" means no result within a fixed 30-minute,
  single-machine compute budget — not proof that a region is impossible to
  index or map.** Two regions encountered during earlier development of
  this methodology (`cand_chr21_347`, `cand_chr21_061`) initially looked
  like failures under a 10-minute budget and turned out to complete
  successfully once given 30 minutes (`cand_chr21_347` is in fact one of
  this analysis's 50 regions — the single highest-complexity region tested,
  complexity_score 2.889 — and completed successfully here in ~26 minutes
  at 99.78% accuracy). So the 16% high-tier failure rate should be read as
  an upper bound under this project's specific resource budget, not a
  proven ceiling; a larger compute budget could shrink it.
- **The exact percentages carry real statistical uncertainty.** The 16%
  figure comes from only 4 failure events out of 25 high-complexity
  regions, and the hazard-vs-failure comparison (4% vs. 12%) from even
  fewer events per cell. The *pattern* is clear — high complexity tier
  fails more often than low; hazard status doesn't predict failure the way
  complexity does — but with this few events, the precise percentages
  should be read as approximate, not as precise population parameters.

### How this design was arrived at

Earlier exploratory work in this project — 4 hand-picked regions spanning
the complexity range (lowest / lower-middle / upper-middle / highest,
short-handed `A`-`D`) — surfaced two problems that directly shaped the
50-region stratified design used here. Full narrative in
`runs/logs/D_highest_substitution_saga.md` and
`runs/logs/haplotype_representativeness_correction.md`;
condensed here:

- **Picking the single highest-complexity window directly was not viable.**
  The top-scoring window in the original 345-window pool,
  `cand_chr21_010` (complexity_score 3.452, over 3x any other window),
  turned out to be **ribosomal DNA** (the 45S rDNA repeat array on chr21's
  acrocentric short arm, containing `RNA5-8SN3` and `MIR6724`) — one of the
  most repetitive sequence classes in the human genome, and its Giraffe
  index never completed even after 30 combined minutes with a stable,
  non-growing memory footprint. Four more attempts at high-complexity
  substitutes (`cand_chr21_347`, `_117`, `_354`, `_329`) also failed at
  various pipeline stages before a working region was found — a direct
  early hint that pipeline tractability, not just accuracy, needed to be
  measured systematically rather than assumed. (`cand_chr21_010` itself
  also turned out to sit inside the pericentromeric corridor later excluded
  by the haplotype-representativeness filter below — consistent with rDNA
  being poorly resolved in HPRC draft assemblies.)
- **A degenerate-window filter looked adequate but validated the wrong
  proxy metric, twice.** First, a z-score corruption bug (see Methods) let
  a degenerate window score as a spurious highest-complexity outlier before
  a pre-filter caught it. Second, and more seriously, the original
  lowest-complexity pick, `cand_chr21_017`, had only 19 of 466 haplotypes
  actually represented — invisible to the existing filter, which counts
  fragmented path pieces rather than distinct haplotypes — and was only
  caught when the user visually inspected its GFA in Bandage. This directly
  motivated the systematic haplotype-representativeness filter in Methods,
  and, once fixed, revealed that the *corrected* lowest-complexity region,
  `cand_chr21_050`, sits on a genuine segmental-duplication hotspot and
  mapped worse than every other region tested at the time — despite scoring
  as the simplest by a wide margin.

Together, these findings motivated two changes that define the current
methodology: (1) building the known-hazard annotation *in advance*, from
public UCSC tracks, instead of discovering problem sequence classes by
accident after each failure; and (2) replacing hand-picked extreme regions
with the properly stratified, 50-region sample this README reports as
primary. Two additional region-specific findings from this phase remain
useful color: `cand_chr21_354`
contains 3 NUMTs (nuclear-mitochondrial pseudogenes) plus a 2,254 bp tandem
array; `cand_chr21_117`'s driver (found only by localizing a node-degree
spike within the window, not from whole-window annotation) is a dense
cluster of 29 overlapping AT-rich microsatellites in a ~2.4 kb stretch.

## Failed analyses

Failures kept and documented rather than silently discarded:

- 22 of 367 candidate windows were degenerate (see "Candidate window
  screening"); excluded before scoring rather than left in to silently
  corrupt the z-score composite via the `NaN`-sigma bug (see "Methods").
- 19 of the remaining 345 windows had severely reduced haplotype
  representation (a contiguous pericentromeric-flanking corridor); excluded
  by the haplotype-representativeness filter (see "Methods" and "How this
  design was arrived at").
- **4 of the 50 final regions** (`cand_chr21_000`, `_046`, `_244`, `_313`)
  did not complete within the 30-minute compute budget — see "Results" for
  the numbers and the resource-budget caveat, and "Technical limitations"
  for the underlying `vg giraffe`/`vg surject` pathology.
- During earlier development of this methodology, 5 consecutive attempts at
  a hand-picked highest-complexity region failed before one succeeded, and
  an earlier haplotype-representativeness bug let an unrepresentative
  window through undetected — see "How this design was arrived at" above
  for the condensed story and the full log files for every attempt.
- `odgi extract -E` failed in every tested configuration (see "Technical
  limitations") — not used in the final pipeline, but the failure itself
  (including one process crash) is a genuine, documented finding about the
  limits of this specific tool/graph combination.
- The `vg` build itself failed 5 times before succeeding (see
  `logs/vg_build_summary.txt`) — every failure and fix is logged there.

## Technical limitations

- **No native macOS `vg` distribution exists.** `vg` is not published for
  osx-arm64/osx-64 via bioconda, is not in Homebrew, and ships only Linux
  binaries/Docker images upstream. It had to be built from source locally
  under `software/vg/` (see `logs/vg_build_summary.txt` for the exact
  toolchain issues hit and how each was resolved — all fixes were additive
  conda-forge packages or environment-variable/build-flag adjustments; no
  system files were modified and no compute was moved off this Mac).
- **Chromosome-level rather than a pre-identified small partition GFA was
  used as the starting graph.** The HPRC release2 PGGB whole-genome build
  also publishes ~1000+ small "community partition" GFAs, which are smaller
  individually, but neither the S3 listing nor the resource README provide a
  partition-to-GRCh38-coordinate mapping, and inspecting on the order of a
  thousand files to locate one covering an arbitrary candidate region was not
  practical within the project timeline. The smallest whole autosome (chr21,
  ~2.4 GB compressed / ~16 GB decompressed GFA) was used instead as a
  tractable stand-in for "partition-level."
- **`odgi extract`'s `-E`/`--full-range` flag was impractical on this graph
  at any tested locus, not only near the centromere.** `-E` is documented to
  need care on "very complex graphs." A near-centromere candidate expanded
  to a ~19 GB `.og` (vs. ~40 MB without `-E`) for a requested 100 kb
  window — nearly the whole 23.7 GB chromosome graph; the same expansion
  (~19 GB `.og` / ~10.6 GB GFA) occurred at every mid-chromosome region
  tested too, and one run crashed outright with a `Bus error` under
  sustained near-total memory pressure. The most likely root cause: `-E`'s
  "sorted order of the graph" span depends on node IDs correlating with
  genomic coordinate order, and `chr21.og` was built without
  `odgi build -s` (topological sort), so the ID order does not stay
  well-localized to each 100 kb window. **All region extraction in this
  project (screening and final regions) was done without `-E`.** Screening
  also disables default subpath-merging (`-d 0`) for speed; final region
  extraction uses the default merging distance, confirmed via `odgi viz`
  inspection to correctly stitch each haplotype's local traversal into ~1
  path per haplotype rather than the tens of thousands of tiny fragments
  `-d 0` produces. Consequence: a haplotype path might in principle be very
  slightly truncated at the window boundary rather than fully "laced" the
  way `-E` would (impractically) attempt — not visibly detectable when
  inspected in practice.
- **Default subpath merging can itself balloon a window far past its
  requested size when it sits adjacent to a large collapsed-repeat
  structure.** Observed independently in at least 4 different regions
  across this project, from 2.5x up to 46x the requested 100 kb. Resolved
  with a `merge_distance=0` override
  (`scripts/shell/extract_region_final.sh`'s 8th parameter), trading a small
  amount of haplotype-path fragmentation for a correctly-sized window. Given
  this recurs unpredictably, `scripts/shell/run_full_pipeline_final_n50.sh`
  now detects and retries it automatically for every region rather than
  relying on a human to notice.
- **`vg giraffe`/`vg surject` exhibit a recurring, unresolved pathological
  hang/crash on specific reads in specific locally-extracted subgraphs, not
  cleanly correlated with graph complexity, memory pressure, or coarse
  genomic position.** First seen as a 21-hour unmonitored hang during a
  since-removed 250 kb pilot extension (Watchdog reporting an implausible
  ~44 TB "memory growth" on one thread); recurred repeatedly during this
  project's region-selection and final-sample runs, with growth figures up
  to ~87 TB reported on a 128 GB machine — physically impossible as true
  memory use, indicating a computational/algorithmic pathology (the tool
  spinning on specific input) rather than a resource-exhaustion problem a
  bigger memory budget would fix. This is the direct cause of the 4 final
  n=50 failures reported under Results, and of several failures during
  earlier region-selection work (see "How this design was arrived at"). No
  root cause was identified; this remains an open, unresolved technical
  limitation of using `vg giraffe`/`vg surject` on PGGB-derived local
  subgraphs at this scale.
- **A 10-20 minute watchdog timeout was too aggressive for some
  genuinely-slow-but-not-hung regions**, directly motivating this project's
  final 30-minute budget — see the resource-budget caveat under Results.
  Every pipeline stage in this project is wrapped in an enforced, logged
  wall-clock watchdog (background PID + `kill -9`, since macOS lacks GNU
  `timeout`) rather than relying on remembering to do so each time, after
  one `vg autoindex` invocation was once launched outside this discipline
  and ran unobserved for 52 minutes before being caught.
- This is a small-scale (single chromosome, 50 regions carried through the
  full pipeline) proof-of-concept; see the analysis spec's own interpretive
  constraints (section 70) for what this analysis does and does not support
  concluding. The full 367-window *screening* pass gives broad coverage of
  chr21's complexity distribution, but only 50 regions were carried through
  indexing/mapping/evaluation.

## Biological limitations

- Primary-analysis reads are simulated error-free from the linear GRCh38
  reference interval of each region, not real sequencing data, so this
  measures graph/mapper behavior under idealized conditions, not real
  Illumina sequencing performance.
- High-complexity regions may have multiple biologically valid mapping
  locations (e.g. segmental duplications); a mismatch to the single simulated
  truth coordinate is not automatically a mapper error. A full repeat-aware
  truth model was out of scope for this project.
- This is not merely a hypothetical concern: the single most extreme
  complexity outlier found during this project's region-selection work,
  `cand_chr21_010` (complexity_score 3.452), was identified via Ensembl
  gene annotation as **ribosomal DNA (rDNA)** — the 45S rDNA repeat array
  on chr21's acrocentric short arm. rDNA is a genuinely, biologically
  repetitive sequence class, not an artifact of this project's complexity
  metric or of PGGB's graph construction — this project's failure to index
  it is a real reflection of a real, well-known hard region of the genome,
  which independently corroborates that the complexity score is at least
  detecting *something* biologically meaningful at its extreme, even though
  (per Results) it does not predict mapping-accuracy degradation the way
  the hazard annotation does.

## Conclusions

Answering the spec's Q1-Q11 (the only window size run in this project was
100 kb — see How to reproduce):

- **Q1 (100 kb index construction):** Mostly yes, but not universally: 46 of
  50 stratified regions completed the full indexing/mapping pipeline; 4
  (16% of the high-complexity tier, 0% of the low tier) did not complete
  within a 30-minute compute budget (see Results for the resource-budget
  caveat on that number). Construction/mapping succeeding is not a safe
  assumption once sampling moves beyond a handful of hand-picked regions.
- **Q2 (250 kb extension):** Not run in this project — out of scope (see
  How to reproduce).
- **Q3 (complexity level where resources increase):** There is no single
  complexity *value* above which failure reliably begins — failures and
  successes both occur throughout the upper half of the range — but
  splitting the population at its **median** complexity_score cleanly
  separates failure risk (0% below the median, 16% above it, subject to the
  statistical-uncertainty caveat in Results). A coarse "top half vs. bottom
  half" threshold is meaningful even though no precise cutoff value is.
- **Q4 (complexity vs. indexing time):** No clean monotonic relationship
  across the full n=50 sample; indexing time is dominated by whether a
  region completes at all (see Q1) rather than by a smooth complexity
  gradient among the regions that do complete.
- **Q5 (complexity vs. peak memory):** Similarly, no clean monotonic
  relationship among completed regions; the more informative split is
  completion vs. non-completion by complexity tier (Q1), not peak memory
  among successes.
- **Q6 (complexity vs. mapping accuracy):** No — this is one of the
  project's more important findings. The single worst-accuracy region in
  the entire 50-region sample, `cand_chr21_050`, is also the single
  *lowest*-scoring region on complexity — it sits on a genuine segmental-
  duplication hotspot invisible to a topology-only score. Mapping accuracy
  is instead well predicted by the independent hazard annotation (see Q9).
- **Q7 (complexity vs. MAPQ):** Same pattern as Q6 — not predicted by
  complexity score; see the hazard-annotation result in Q9 instead.
- **Q8 (size vs. topology):** this project held window size fixed (100 kb)
  across all regions by design, so it speaks only to the topology side.
  Separating size and topology fully would require the size comparison in
  spec section 58, which was not run here.
- **Q9 (most promising predictor metric):** Answered rigorously via the
  stratified n=50 design — and the answer is more precise than
  "no metric works": **`complexity_score`
  and the known-hazard annotation each predict a *different* outcome, and
  neither predicts the other's outcome well.** Splitting `complexity_score`
  at its population median cleanly separates outright pipeline failure (0%
  below median vs. 16% above, n=25 each) — so the composite score *does*
  have real predictive value for tractability, just not a precise threshold
  value (see Q3), and just an upper-bound-under-our-budget value at that
  (see the resource-budget caveat in Results). Separately, `hazard_flag`
  (built from two established UCSC sequence-annotation tracks, calibrated
  to flag exactly the loci this project independently discovered by hand)
  cleanly separates mapping-accuracy degradation among regions that *do*
  complete (minimum accuracy 99.53% without the flag vs. 98.69% with it),
  while `complexity_score` does not. This makes intuitive sense in
  hindsight: pipeline tractability plausibly tracks raw local graph
  scale/tangle (what the composite score measures, however roughly), while
  mapping accuracy degrades from genuine sequence ambiguity (multiple
  near-identical loci a read could originate from) — a property of sequence
  *content*, invisible to a metric computed purely from one window's local
  topology.
- **Q10 (biggest technical limitation):** Three, in practice tied. (1)
  Getting `vg` to build at all on macOS (5 failed attempts before success)
  and discovering that `odgi extract -E` was impractical on this graph. (2)
  The recurring, unresolved `vg giraffe`/`vg surject` pathological
  hang/crash (implausible multi-TB "memory growth," specific reads,
  specific local subgraphs) responsible for all 4 final-sample failures and
  several failures during earlier region-selection work — compounded by the
  fact that at least 2 comparable-looking "failures" turned out, on retest
  with a longer timeout, to be genuinely slow rather than hung, meaning the
  true failure rate under an even more generous budget is unknown and could
  be lower than 16%. (3) A degenerate-window filter that looked adequate but
  validated the wrong proxy metric twice: first a z-score `-inf`/`NaN`
  corruption bug, then a haplotype-representativeness bug that let a
  19-of-466-haplotype window through undetected until a human visually
  inspected the GFA in Bandage. Of the three, (3) is arguably the most
  methodologically significant, since it silently corrupted a headline
  result rather than causing a visible crash.
- **Q11 (worth extending genome-wide?):** Yes. Three concrete
  recommendations for a genome-wide extension, all directly informed by
  this project's own experience: (1) **budget for outright pipeline
  failure as its own outcome, with a resource-aware definition of
  "failure"** — this project's own n=50 sample had a 16% failure rate among
  above-median-complexity windows under a 30-minute budget (0% below
  median), and that figure is itself an upper bound, not a hard ceiling; (2)
  **screen candidate windows against known-hazard sequence annotation
  (segmental duplications, large tandem repeats, rDNA) before trusting a
  local-topology complexity score's ranking** — this project's own
  `hazard_flag`, built from two public UCSC tracks, predicted
  mapping-accuracy degradation far better than `complexity_score` did, and
  is cheap to compute in advance; (3) **use a generous (>=30 minute) per-
  region timeout, and treat the resulting failure rate as budget-dependent**
  — this project caught 2 false failures caused purely by an insufficiently
  patient watchdog, and cannot rule out that some fraction of its own
  reported 16% would resolve given more time or compute.

Scientific-interpretation framing: the 100 kb results sit closest to
*"mapping remained mostly robust across the tested regions once indexes
could be constructed, but 'robust,' 'once constructed,' and the exact
failure rate itself are each carrying a specific caveat this project can
now name precisely."* First, mapping-accuracy degradation among regions
that *do* complete is **not** explained by complexity score at all (the
single worst-accuracy region of all 50 tested is also the single
*lowest*-scoring one) — it is explained by the independently-checkable
known-hazard annotation instead. Second, "once indexes could be
constructed" is doing real work: above-median-complexity windows failed
outright 16% of the time under this project's compute budget (vs. 0% below
the median) — a real, replicated risk, not a one-off. Third, and the
caveat that qualifies the second: that 16% is a property of a fixed
30-minute, single-machine budget, not a proven mathematical ceiling, and
rests on only 4 failure events — a larger compute budget or a longer
timeout might resolve some fraction of it, and the precise percentage
should not be over-read. The project's most defensible, final conclusion is
therefore not "this complexity score is simply unreliable" but the more
precise: **local-topology complexity score reliably predicts whether the
pipeline completes within a given compute budget, but not mapping accuracy
among completed regions — for that, checking directly for known-hazardous
sequence classes (rDNA, large tandem repeats, NUMTs, segmental
duplications) works better and should be done explicitly rather than
assumed to correlate with the score.**

## How to reproduce

1. `cd PGGB_short_reads_mapping` (this directory).
2. Ensure the `bh26` conda environment exists and is active (see
   `results/bh26_environment.yml` for an exact re-creation spec, or
   `results/bh26_environment_from_history.yml` for the minimal explicit
   package list).
3. Run the scripts under `scripts/shell/` and `scripts/python/` in the order
   described in Methods above, with paths pointed at
   `runs/`:
   - `scripts/python/screen_candidate_windows.py --all` to enumerate every
     window, then `scripts/shell/run_screen_batch_safe.sh` to screen them
     all with an enforced per-window timeout.
   - Apply the degenerate-window and haplotype-representativeness filters
     (see "Candidate window screening"), then
     `scripts/python/select_representative_regions.py` to compute
     `complexity_score` on the resulting 326 windows.
   - Build the known-hazard annotation: download UCSC `genomicSuperDups` and
     `simpleRepeat` for the whole chromosome and intersect locally against
     `results/graph_stats.tsv`'s 326 valid windows (see
     `runs/logs/hazard_stratified_robustness_analysis.md`
     for the exact thresholds).
   - Rank-sample 50 regions within each of the 4
     `{complexity_tier} x {hazard_flag}` cells into a manifest
     (`region_id, chrom, start, end`), then run
     `scripts/shell/run_full_pipeline_final_n50.sh <manifest.tsv> <run_dir>
     12 1800` (a **30-minute**, not 10-20 minute, per-region timeout — see
     "Technical limitations" for why this matters) to process the whole
     batch unattended, with automatic detection/retry of the default-merge
     extraction blowup pathology and mapping evaluation included end to end.
   - Build the final table:
     `scripts/python/build_indexing_stats.py` then
     `scripts/python/build_combined_results_n50_final.py`, and regenerate
     figures with `scripts/python/make_figure1_complexity_distribution.py`
     and `scripts/python/make_figure2_hazard_vs_complexity.py`.
   - Long-running steps are designed to be launched inside `tmux` with an
     external wall-clock watchdog (see each script's header comment and
     "Technical limitations") so they survive independently of any
     particular shell/agent session and cannot hang unmonitored.
4. All raw inputs are re-downloadable from the URLs in
   `metadata/input_files.tsv`; nothing needed to reproduce this analysis lives
   outside this project directory or the `bh26` conda environment.

**Scope note:** this project ran the 100 kb primary analysis only. The
optional 250 kb same-locus size comparison and the 500 kb /
Minigraph-Cactus extensions (both explicitly optional in the analysis spec)
were not run.
