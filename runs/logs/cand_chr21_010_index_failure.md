# cand_chr21_010 (extreme-complexity 100 kb window) — Giraffe index infeasible

## Summary
`vg autoindex --workflow sr-giraffe` could not complete for this region within
a combined 30 minutes across two attempts (10 min, then 20 min). This is
recorded as a genuine "index construction infeasible at extreme local
complexity" finding, not retried further.

## The region
- chr21:8,334,199-8,434,199 (100 kb, 0-based start), immediately adjacent to
  a cluster of degenerate/collapsed-node windows spanning chr21:7.9-8.6 Mb
  and chr21:10.9-12.8 Mb (see results/graph_stats_excluded_degenerate.tsv).
- Extracted (without subpath merging, `-d 0`, to keep the correct 100 kb
  window size -- default merging ballooned this same region to 4.6 Mb / 46x
  the requested size on a first attempt, since it chased fragments into the
  adjacent collapsed structure).
- Raw graph stats at the correct 100 kb size: 24,108 nodes, 30,394 edges,
  **2,719,131 paths**, **205,571,872 steps**. For comparison, the pilot's
  most complex 100 kb region had ~78,000 paths and ~1.7M steps -- this
  region has ~35x more paths and ~120x more steps at the *same* window size.

## What happened during indexing
- GBWT construction (path/haplotype indexing) and the Giraffe distance index
  both completed quickly and normally (~37s and a few seconds respectively;
  `cand_chr21_010.dist`, 12.7 MB, was produced both times).
- The **minimizer/zipcode index** (`build_minimizer_index`, k=29 w=11) never
  completed. Both attempts stopped producing any further log output after
  "Cached payloads in 0.0158889 seconds" -- identical output on both runs,
  suggesting whatever comes next (the actual minimizer enumeration across
  2.7M path traversals) does not print incremental progress and simply did
  not finish in the allotted time.
- **Memory was stable and low throughout both attempts** (~1-2 GB tracked
  RSS, no growth, no swap) -- this is explicitly *not* the same failure mode
  as the earlier 250 kb `vg surject` incident (which showed runaway,
  implausible memory growth). This is instead simple computational
  intractability: the minimizer-indexing algorithm's cost apparently scales
  with total path-traversal length (205M steps here) in a way that makes it
  impractically slow at this scale, not a bug or a resource leak.

## Biological identity: ribosomal DNA (rDNA)

Looking up GRCh38 chr21:8,334,199-8,434,199 against Ensembl
(`overlap/region` REST endpoint, `feature=gene`) identifies this window as
part of the **45S ribosomal DNA (rDNA) repeat array** on chr21's acrocentric
short arm (the nucleolar organizer region, NOR): it contains `RNA5-8SN3` (a
5.8S ribosomal RNA gene) and two `MIR6724` miRNA genes, both known to lie
within the rDNA repeat unit's spacer sequence, plus several lncRNA/TEC/
misc_RNA annotations and no protein-coding genes. rDNA is one of the most
repetitive, poorly-resolved sequence classes in the human genome (GRCh38
itself only partially/fragmentarily represents it), so a PGGB graph
representing near-identical tandem repeat copies across 466 haplotypes as a
densely tangled bubble structure here is entirely consistent with known
biology, not an artifact of this project's methods. This explains why this
single window was such an extreme complexity_score outlier (3.452, next
highest was ~0.9) -- but note it does *not* explain the other 4 substitute
failures encountered later in the `D_highest` search (`cand_chr21_347`,
`_117`, `_354`, `_329`), none of which are rDNA/NOR and none of which scored
anywhere near this extreme (see
`runs/logs/D_highest_substitution_saga.md`).

## Disposition
- Not used as the working "D_highest" region for the main 4-region
  comparison. Five substitute candidates were attempted in total
  (`cand_chr21_347`, `_117`, `_354`, `_329`, in that order after this one),
  and all five also failed at some pipeline stage -- see
  `runs/logs/D_highest_substitution_saga.md` for the
  full account. `cand_chr21_330` (the pilot's own previously-validated
  `D_highest` region) was ultimately adopted instead.
- `cand_chr21_010`'s extraction *was* completed successfully (see
  `work/extracted_final/cand_chr21_010/odgi_stats.tsv`) and is kept as a
  standalone data point: proof that, among the 367 windows screened,
  extraction/graph-construction remained possible even for this extreme
  case, but downstream Giraffe index construction did not.
- This is treated as a genuine result relevant to the project's central
  research question -- a concrete example of "at what level of local
  complexity does short-read mapping become computationally difficult" --
  not an operational footnote to be hidden.
