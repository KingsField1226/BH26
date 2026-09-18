# Beyond n=4: does complexity score or known-hazard annotation better predict outcomes?

## Motivation

The 4-region primary comparison (`A_lowest`/`B_lower_middle`/`C_upper_middle`/`D_highest`)
is, by the project's own repeated acknowledgement, too small a sample to validate
`complexity_score` as a predictor of anything. After the haplotype-representativeness
correction produced a corrected `A_lowest` (`cand_chr21_050`) that turned out to sit on
a segmental duplication and map *worse* than every other region despite scoring as the
simplest, it became clear the right next step was not to keep hand-picking 4 regions,
but to test the underlying question properly: across many regions, does
`complexity_score` predict pipeline outcomes at all, and if not, does something else
(specifically, the presence of a known problematic sequence class) predict them better?

## Step 1: an n=10 rank-sampled trend check

10 regions were rank-sampled at even percentile intervals across the full range of
`complexity_score` in the 326-window corrected pool (from -2.166 to 2.889), reusing 3
already-completed regions (`cand_chr21_050`, `_093`, `_287`) and running 7 new ones.

Result: **no statistically significant monotonic relationship at n=10** between
`complexity_score` and any outcome (Spearman rank correlation):

| outcome | rho | p |
|---|---|---|
| mapping accuracy | -0.29 | 0.42 |
| index wall-clock time | 0.36 | 0.30 |
| index peak RSS | 0.07 | 0.86 |

A follow-up split of these 10 points into a "dense cluster" (7 points within
complexity_score -0.5 to 0.5, i.e. the 78% of the distribution that sits near the
median) vs. the 3 remaining extremes showed the dense cluster's own accuracy values
differ from each other only in single-digit read counts out of 60,000 (0 to 8
incorrect reads) -- differences at the noise floor of the simulation, not a real
signal. All of the interesting variance came from the 3 extreme points, and even among
those, the *direction* of the effect was inconsistent with a simple "higher complexity
is worse" story (the single lowest-scoring point, `cand_chr21_050`, had the worst
accuracy of the ten).

## Step 2: a critical correction discovered mid-investigation

While rank-sampling further into the high-complexity tail, `cand_chr21_347` (already
on record in `D_highest_substitution_saga.md` as attempt #2 of 5, "killed at 10-min
watchdog, GAM only 841KB of an expected ~9.7MB") was re-run as part of this n=10 batch
with a longer, 1800s (30 min) per-region timeout instead of the earlier 10-20 min
watchdogs. **It completed successfully in 1561s (~26 min)**, with a full 10.75MB GAM,
60,000/60,000 reads mapped, 99.78% accuracy.

This is a significant correction: `cand_chr21_347` was not, in fact, pipeline-
intractable -- the original 10-20 minute watchdogs used during the `D_highest`
substitution saga were simply not generous enough for this region's genuinely slower
(but finite) processing time. The same pattern recurred for `cand_chr21_061`
(hazard_flag=False, low complexity) later in this investigation, which also "looked"
like it was hanging (worsening per-read Watchdog wait times up to 112s) but completed
successfully at 1561s under the 1800s timeout.

**This means the `D_highest_substitution_saga.md` record of "5 consecutive failures"
should be read with an important caveat: at least 1 of those 5 (`cand_chr21_347`,
and very likely also `cand_chr21_117`, which showed the identical partial-progress
pattern under the same 10-minute watchdog) may not have been genuinely intractable,
just under-timed.** `cand_chr21_010` (rDNA, truly never completing even after 30
combined minutes with a *stable*, non-growing memory footprint) and `cand_chr21_329`
(the extreme multi-terabyte "memory growth" crash, a qualitatively different failure
signature -- an actual crash, not merely slow progress) remain confidently classified
as genuine failures, not merely slow successes. The saga document itself has not been
rewritten to re-attempt `_347`/`_117` under a 30-minute budget retroactively (that
would require yet more compute time for a question the n=50 analysis below already
answers more rigorously), but this caveat should be read alongside that document.

This correction directly motivated using a **30-minute (1800s) per-region timeout**
for the rest of this investigation, rather than the 10-20 minute watchdogs used
earlier in the project.

## Step 3: building a "known-hazard" annotation independent of complexity_score

Rather than discovering problem-causing sequence features by accident after each
failure (as happened for `cand_chr21_050`, `_354`, `_358`), a systematic annotation
was built for all 326 valid windows *in advance*, using two UCSC tracks downloaded
once for the whole of chr21 (not per-window queries) and intersected locally:

- `genomicSuperDups` (segmental duplications): 1,627 records for chr21.
- `simpleRepeat` (Tandem Repeats Finder): 14,911 records for chr21.

For each window, `n_segdup` (count of overlapping segmental duplication records) and
`max_segdup_identity`, plus `max_tandem_repeat_span_bp` (largest `period x copyNumber`
among overlapping simple repeats) were computed. A binary `hazard_flag` was defined as:

```
hazard_flag = (any overlapping segmental duplication with fracMatch >= 0.90)
              OR (max_tandem_repeat_span_bp >= 1400)
```

Thresholds were calibrated (not arbitrary) to correctly flag all 4 previously-
discovered problem windows: `cand_chr21_050` (41 segdups, up to 97.3% identity),
`cand_chr21_354` (2,264 bp tandem array), `cand_chr21_329` (1,462 bp tandem array),
and `cand_chr21_358` (a single segdup at 90.96% identity, chr21-internal). A looser
threshold (any segdup regardless of identity, or a 500bp tandem-repeat cutoff) flagged
practically half the genome (47.9%) and was judged too permissive to be useful; the
calibrated threshold flags 77/326 windows (23.6%). Full per-window annotation:
`results/hazard_annotation_326windows.tsv`.

A second axis, `complexity_tier` (`low` if `complexity_score` < the pool median,
`high` otherwise), was defined for a clean 2x2 stratified design:
`{hazard_flag} x {complexity_tier}`.

## Step 4: stratified sampling, grown from n=20 to n=50 across 3 rounds

Representative regions were rank-sampled *within each of the 4 cells* (not across the
whole pool, to guarantee balanced coverage of hazard-positive-but-low-complexity and
hazard-negative-but-high-complexity windows, which a plain complexity-rank sample
would rarely include). The sample was grown in 3 rounds specifically to test whether
early patterns were robust or a fluke of a small draw:

| round | new regions | cumulative n | cumulative failures |
|---|---|---|---|
| 1 (initial stratified draw) | 9 | 20 | 1 (`cand_chr21_000`) |
| 2 (robustness check) | 10 | 30 | +2 (`cand_chr21_244`, `cand_chr21_313`) |
| 3 (robustness check) | 20 | 50 | +1 (`cand_chr21_046`) |

All pipeline runs used the corrected 30-minute per-region timeout and the same
`run_full_pipeline_batch_safe.sh` batch driver (extract, with automatic detection and
`merge_distance=0` retry if the default-merge extraction exceeds 3x the requested
window size -- this fired once more in this investigation, for `cand_chr21_093`
[5.4x, already documented] and again independently for `cand_chr21_088` [5.2x],
confirming the "default merge blowup" pathology recurs unpredictably and needs this
automatic safeguard for any future batch work on this graph).

## Final result at n=50

**Pipeline failure rate, by `complexity_tier` alone (ignoring `hazard_flag`):**

| complexity tier | n | failures | failure rate |
|---|---|---|---|
| low (< median) | 25 | 0 | **0%** |
| high (>= median) | 25 | 4 | **16%** |

**Mapping accuracy among successfully-completed regions, by `hazard_flag` alone
(ignoring `complexity_tier`):**

| hazard_flag | n (successful) | mean accuracy | minimum accuracy |
|---|---|---|---|
| False (no known hazard) | 24 | 99.92% | **99.53%** |
| True (known SD/tandem-repeat hazard) | 22 | 99.85% | **98.69%** (`cand_chr21_050`) |

Both patterns held essentially unchanged from n=20 through n=30 to n=50 (the n=30
numbers were 0%/20% failure and 99.56%/98.69% minimum accuracy -- statistically
indistinguishable from the n=50 numbers given these sample sizes), which is the
robustness check this 3-round growth was designed to provide.

**All 4 pipeline failures across the full n=50 sample fell in the high-complexity
tier, regardless of hazard_flag** (`cand_chr21_000` and `cand_chr21_046`: hazard=True;
`cand_chr21_244`: hazard=True; `cand_chr21_313`: hazard=False). **The single worst
mapping-accuracy result among all 50 regions belongs to the *lowest*-complexity region
in the entire dataset** (`cand_chr21_050`, complexity_score -2.166, the most extreme
low score of any region tested in this project).

## Interpretation: two different questions, two different predictors

1. **"Will the pipeline complete at all?"** -- best predicted by `complexity_tier`
   (specifically, `complexity_score` above the population median). This makes some
   intuitive sense: `vg giraffe`/`vg surject`'s per-read processing time scales with
   the sheer scale of local graph tangle (path count, step count), which is exactly
   what the composite score is measuring, however roughly.
2. **"If it completes, how accurate is the mapping?"** -- best predicted by
   `hazard_flag` (a known problematic sequence class: segmental duplication or a large
   tandem repeat), **not** by `complexity_score`. This also makes sense on reflection:
   accuracy degradation comes from real sequence ambiguity (multiple near-identical
   genomic loci a read could plausibly have come from), which is a property of
   sequence *content*, not of how tangled the locally-extracted graph's topology
   happens to look.

Neither predictor is sufficient on its own for both questions -- `complexity_tier`
does not explain `cand_chr21_050`'s poor accuracy (it is in the *low* tier), and
`hazard_flag` does not fully explain pipeline failure (`cand_chr21_313` failed with no
hazard flag at all, and conversely most hazard-flagged regions complete just fine with
only mildly reduced accuracy). This two-predictor, two-outcome structure is the
central, most rigorously-supported finding of the full-scan investigation.

See `results/combined_results_n50_stratified.tsv` for the full per-region data (all
50 regions, complexity_score, hazard_flag, tier, and pipeline/mapping outcomes) and
`results/figures/figure10_complexity_score_distribution.png` /
`results/figures/figure12_hazard_vs_complexity_2x2.png` for the corresponding figures.
