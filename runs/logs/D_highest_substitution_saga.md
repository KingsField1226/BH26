# D_highest region selection: five failed candidates before settling on cand_chr21_330

## Summary

For the full 367-window scan, after screening produced 345 valid (non-degenerate)
candidate windows and scoring them by `complexity_score`, the intent was to pick the
single highest-scoring window as `D_highest`, exactly as was done for the 10-window
pilot. In practice, **five successive high/moderate-complexity candidates all failed**
somewhere in the indexing/mapping/surject pipeline, each with the same pathological
signature previously seen in the (now-removed) 250kb pilot extension: `vg` reporting
absurd, physically-impossible per-thread "memory growth" figures (gigabytes to tens of
terabytes on a 128GB machine) while stalling on a small number of specific reads, well
before any real memory pressure develops. This is documented here as a first-class
finding, not hidden as an inconvenience: **the current 4-metric complexity score does
not reliably predict vg Giraffe/surject computational tractability**, and the failure
mode appears sensitive to genomic position at a finer resolution than 100kb windowing
resolves.

**Attempt 1 has a clear-cut biological cause; attempts 2-5 all show a
plausible repeat-driven contributing factor, though attempt 3's was only
found after localizing it precisely (see below) rather than from a
whole-window annotation pass.**
`cand_chr21_010` (chr21:8,334,199-8,434,199) is ribosomal DNA (rDNA) -- the
45S rDNA repeat array / nucleolar organizer region on chr21's acrocentric
short arm, identified via Ensembl gene annotation (`RNA5-8SN3`, `MIR6724-3`,
`MIR6724-4`; see `logs/cand_chr21_010_index_failure.md` for detail). rDNA is
a genuinely, biologically extreme repeat class, which is a satisfying
explanation for why this window's complexity_score (3.452) was such an
outlier.

A follow-up annotation check (Ensembl genes + UCSC `genomicSuperDups` +
UCSC `simpleRepeat`/TRF) on the other 4 attempts found: **no segmental
duplications overlap any of the 4 regions**, and all 4 turned out to contain
substantial tandem-repeat (microsatellite/minisatellite) arrays that plausibly
contribute to graph tangling, though none are as extreme as rDNA and none of
this is proof that these specific repeats are what `vg` stalls on:

- `cand_chr21_347` contains the **CSTB** gene (progressive myoclonic epilepsy
  type 1 is caused by a pathogenic dodecamer-repeat expansion in CSTB's
  promoter) plus several TRF-called tandem arrays up to 47.8 copies of a
  20 bp unit (~980 bp span) and 13.3 copies of a 40 bp unit (~540 bp span).
- `cand_chr21_354` contains **3 NUMTs** (nuclear-mitochondrial pseudogenes:
  `MTND6P21`, `MTCYBP21`, `MTND5P1` -- insertions of mitochondrial sequence
  into the nuclear genome, known to vary in copy number/presence across
  individuals) plus the largest tandem array found in any of the 4 regions:
  **566 copies of a 4 bp unit spanning 2,254 bp**, and an additional
  213 bp-period array (8 copies, 1,729 bp span).
- `cand_chr21_329` contains ordinary protein-coding genes (`PRDM15`,
  `RIPK4`) but has two large, overlapping tandem arrays (174 copies of a
  6 bp unit and 73 copies of the same region's 15 bp higher-order unit,
  together spanning chr21:41,805,979-41,807,094, ~1.1 kb) plus a 31-copy
  47 bp array (chr21:41,813,680-41,815,153, ~1.5 kb) -- both located in the
  back third of the window, i.e. toward the boundary shared with the
  successful `cand_chr21_330` (see "The striking part: adjacency" below).
- `cand_chr21_117` is a **gene desert** (one processed pseudogene,
  `ZNF799P`), and an initial whole-window annotation pass found nothing more
  than an "ordinary" 76-copy TA dinucleotide microsatellite -- but this was
  a false negative from evaluating the whole window's annotation in
  aggregate, not an absence of a driver. **Update, after visually inspecting
  the extracted GFA in Bandage** (the user requested the raw GFA files for
  all 4 regions specifically to do this): the graph has one sharply
  localized tangle, not a diffuse pattern. Per-base node-degree profiling
  along the GRCh38 reference path (`odgi degree -i cand_chr21_117.og -R
  <path-list> -D`, then scanning the per-base degree vector) pinpointed it
  precisely: degree jumps from a baseline of ~4 to as high as **76** across
  a narrow ~3 kb stretch at local offset ~52,000-55,000 bp within the
  window, i.e. absolute **chr21:20,586,138-20,589,147**. UCSC's
  `simpleRepeat`/TRF track shows this ~2.4 kb interval is not one repeat but
  a **dense cluster of 29 overlapping/nested AT-rich short tandem repeats**
  (mostly TA/AC/TG dinucleotides; copy numbers up to 76, 53, 22, 21.5, 18,
  17, 16, 15, plus several 16-44 bp compound units) -- a classic unstable
  microsatellite cluster, exactly the kind of per-haplotype repeat-length
  variation that produces many alternate paths through a PGGB bubble.

**Lesson learned from this correction:** a whole-window annotation summary
can miss a driver that is real but spatially concentrated; per-base
node-degree profiling (`odgi degree -D`) plus visual Bandage inspection is a
more reliable way to localize *where* in a window the graph's complexity
actually lives before querying annotation tracks at that specific position.

None of this is a mechanistic proof that `vg giraffe`/`vg surject` actually
stall *on* these specific repeat loci -- confirming that would require
isolating and profiling the exact offending read(s) inside `vg`, which was
judged out of scope (it would mean re-running more of the same
multi-terabyte-hang-prone jobs this project has already decided to stop
chasing). This is an annotation-correlation check, not a root-cause fix.

## Attempts, in order

| # | region_id | position (chr21) | complexity_score | percentile (of 345) | Stage failed | Signature | Log(s) |
|---|---|---|---|---|---|---|---|
| 1 | cand_chr21_010 (rDNA / NOR -- see above) | 8,334,199-8,434,199 | 3.452 (max) | 100% | `vg autoindex` minimizer/zipcode step | Two attempts (10-min, 20-min watchdog), memory stable ~1-2GB both times -- not a leak, the step itself simply never completes at this scale (2.7M paths / 205M steps) | `logs/cand_chr21_010_index_failure.md`, `logs/index_d10.log`, `logs/index_d10b.log` |
| 2 | cand_chr21_347 | 43,734,199-43,834,199 | (moderate) | -- | `vg giraffe` mapping itself | Killed after ~12min wall / 75min CPU; retried with explicit 10-min watchdog, still timed out; GAM only 841KB of an expected ~9.7MB | `logs/map_347.log`, `logs/map_347_watchdog.log`, `logs/map_347_memory.log` |
| 3 | cand_chr21_117 | 20,534,199-20,634,199 | (moderate) | -- | `vg giraffe` mapping itself | Same pattern; 10-min timeout; GAM 7.6MB (more progress, still incomplete) | `logs/map_117.log`, `logs/map_117_watchdog.log` |
| 4 | cand_chr21_354 | 44,434,199-44,534,199 | (moderate) | -- | `odgi extract` (final, default-merge) | Extraction ballooned to ~2.5x the requested 100kb window (adjacent large collapsed-repeat structure pulled in by default 300kb subpath merging); the follow-on `vg autoindex` was then launched in a plain foreground command (not tmux/watchdog-protected) and ran unmonitored for 52 minutes before being discovered and force-killed. No persisted log survives for this specific incident -- it predates `run_full_pipeline_100kb.sh` and the mandatory tmux+watchdog convention that this incident itself motivated adopting project-wide from this point on. | (none persisted; see README "Technical limitations") |
| 5 | cand_chr21_329 | 41,734,199-41,834,199 | 0.640 | ~92% | `vg surject` | Extraction, indexing, read simulation, and `vg giraffe` (exit 0, 10.75MB GAM) all completed in 255s. `vg surject` then logged Watchdog warnings with the most extreme growth figures seen in this project: "Thread 1 ... 66593341440 kb memory growth" (~66.6 TB) and "Thread 7 ... 86718693376 kb memory growth" (~86.7 TB). The process crashed (consistent with a real allocation failure on an actual 128GB machine following an attempted ~TB-scale allocation), leaving a 171-byte (empty) `.surjected.bam`. tmux session ended on its own; no orphaned process remained. | `logs/pipeline_329.log`, `logs/pipeline_329_watchdog.log`, `logs/pipeline_329_memory.log` |

## The striking part: adjacency

`cand_chr21_329` (41,734,199-41,834,199) sits **immediately adjacent** to the pilot's
own `D_highest` region, `cand_chr21_09` / `cand_chr21_330` (41,834,199-41,934,199) --
a 100kb window sharing a boundary. The pilot's region has been run successfully twice
now (once in the original 10-window pilot, once again in this full scan as
`cand_chr21_330`), producing byte-identical extracted-graph statistics both times
(graph_bp=102719, node_count=5378, path_count=466, accuracy=0.99558). Its immediate
neighbor, at a *lower* complexity_score (0.640 vs 0.445) and comparable raw stats,
crashes `vg surject` catastrophically. Complexity score and even coarse genomic
position are therefore not sufficient predictors of tractability here; whatever
drives the pathology is more fine-grained than 100 kb window placement.

A follow-up check found no segmental duplication overlapping `cand_chr21_329`,
but it does contain two large tandem-repeat arrays (a ~1.1 kb array of a
6/15 bp repeat unit and a ~1.5 kb array of a 47 bp unit; see above) sitting in
the back third of the window, i.e. closest to the shared boundary with
`cand_chr21_330`. This is a plausible, though unconfirmed, candidate for why
the two windows behave so differently: depending on exactly where the 100 kb
cut falls relative to a repeat array, the same underlying repeat structure can
end up almost entirely inside one window's graph or split across a boundary in
a way that changes the extracted subgraph's topology substantially.

## Resolution

After 5 consecutive failures across complexity scores from 0.445 to 3.452 and two
different failure stages (indexing, mapping, surject), the search for a higher-scoring
`D_highest` was stopped. `cand_chr21_330` (the pilot's own proven-working coordinates,
complexity_score=0.445, 87.5th percentile of the 345-window pool) was adopted as
`D_highest` for the full-scan 4-region dataset. This is a deliberate choice to end on
a *validated working* highest-complexity-tier region rather than continue an
open-ended search, and the failed attempts are kept as a primary finding of this
analysis (see README "Technical limitations" and "Conclusions").
