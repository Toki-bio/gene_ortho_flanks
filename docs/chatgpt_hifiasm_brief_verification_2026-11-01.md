# Verification of the ChatGPT hifiasm brief (2026-11-01)

Brief: `chatgpt_hifiasm_brief_2026-11-01.md` (saved verbatim). Source: hifiasm 0.25.0, commit `ec9a8b2`, clone read locally. Assembly: our trizona hifiasm 0.25 Hi-C run, therioserver:/home/toki/gof/asm25/ (logs `hifi.log`, `hic.log`; graphs `tri.hic.p_utg.gfa`, `tri.hic.hap{1,2}.p_ctg.noseq.gfa`).

## Q6. Does our run use the mmhap route (`gen_mmhap_t` / `infer_mmhap_copy`)? No. Confirmed.

- `Overlaps.cpp:39640-39649` chooses the Hi-C output path by ploidy:
  `if (asm_opt.polyploidy <= 2) output_hic_graph(...); else output_hic_graph_mmhap(...);`
  Our run is diploid (default ploidy), so `output_hic_graph` runs.
- `output_hic_graph` calls `hic_analysis(..., NULL /*rh*/, ...)` (`Overlaps.cpp:16316`); with `rh == NULL`, `hic_analysis` calls `hic_short_align`, not `hic_short_align_mmhap` (`hic.cpp:17717-17718`). Only `output_hic_graph_mmhap` passes `&rh` (`Overlaps.cpp:16577`).
- `gen_mmhap_t` always prints `[M::stat] # k-copy bases` (`hic.cpp:17464`). Neither of our logs contains such a line.
- So `infer_mmhap_copy` (`Overlaps.cpp:15688`) and `purge_phase` were not executed for this assembly; the brief's hypothesis 3.2-3.3 describes the polyploid (`--n-hap` > 2) path.

## What did run (Confirmed from the log)

`hic.log`: `[M::purge_dups] homozygous read coverage threshold: 25` and `[M::adjust_utg_by_trio] primary contig coverage range: [21, infinity]`. The latter is printed in `Overlaps.cpp:21006-21022` (lower bound = 0.85 x 25), on the path where `if_primary_unitig` (`Overlaps.cpp:18895-18962`) protects a unitig from being removed from a haplotype (`Overlaps.cpp:20228-20240`).

## New genome-wide test on hifiasm's own output (Confirmed)

`hifi/utg_rule_test.py`: every p_utg unitig >= 20 kb; rd tag, GC, and whether its reads ended in both final haplotypes (A lines). Output therioserver:/home/toki/gof/ends/utg_rule_test.tsv.

| rd | share of length placed in both haplotypes |
|---|---|
| <= 15 | 0-2 % |
| 16 | 25 % |
| 17 | 74 % |
| 18-20 | 82-93 % |
| >= 21 | 95-100 % |

Placement in both haplotypes follows a sharp coverage switch at rd ~16-17. Near the switch GC matters at equal rd: GC >= 48 % unitigs are placed in both 72 % (rd 14) and 62 % (rd 15) of the time, GC < 42 % unitigs 0 %. So hifiasm's internal measure (per-read exact-overlap coverage in `if_primary_unitig`) is not the rd tag itself; the category-1 pieces (GC 49-51 %, rd 13-15, e.g. R01 utg004194l rd 13, R07 utg005092l rd 15, R04 utg007335l rd 12) are GC-rich two-copy unitigs that fell below it.

## Answers to the brief's questions

1. `infer_mmhap_copy` is not used here (Q6). The relevant measure is per-read exact-overlap coverage in `if_primary_unitig`. Confirmed (code), GC depression of coverage measured (k-mer and read depth).
2. Yes: two-copy unitigs at rd 13-15 end up in one haplotype; genome-wide, unitigs at rd <= 15 are almost never placed in both. Confirmed (table).
3. Threshold: 0.85 x H = 21 on per-read exact-overlap coverage, >= 80 % of reads or of overlap bases (`Overlaps.cpp:18895-18962`, `21006-21022`); observed switch at rd 16-17. Confirmed (code + data); exact per-read values for each case not recomputed.
4. A unitig not protected as primary is removed from the haplotype whose reads it lacks (`Overlaps.cpp:20228-20240`). Confirmed (code).
5. Hi-C does not restore a removed unitig in this path; no code re-adds an unprotected unitig because of balanced Hi-C signal (code reading; strongly supported).
6. No (above). Confirmed.
7. Before scaffolding: the second copy is missing from `tri.hic.hap{1,2}.p_ctg` already; the stretch is one unitig placed in one haplotype. Confirmed for the 10 traced regions (graph_trace.py, utg_by_bin.py).
8. The N is written by Hi-C scaffolding (YaHS in our run; the deposited assembly lists "hifiasm, YaSH, PretextView"), across a hole that hifiasm left. Confirmed for our run; for the deposited run inferred from identical behaviour.
9. Falsification evidence would be reads joining the flanks directly (copy B lacking the piece) or a parallel copy-B unitig; none found at 10 regions and 40 sampled gaps. Independent two-copy evidence: no direct-join or other-sequence reads, flanks present in both haplotypes, betulina has the piece. Copy number itself rests on GC-calibrated k-mers.
10. Alternative for the 340 "unreadable" gaps: one copy's coverage fades (allele-specific HiFi dropout), a different mechanism, kept separate.

## Errors in the brief

- Its central mechanism (mmhap copy-number inference) is not on the diploid execution path.
- Its gap table drops no category but mis-sums: 76 + 340 + 65 + 47 + 72 + 36 + 29 + 10 + 56 = 731; "the first eight sum to 666" is wrong (they sum to 666 only without the 65; the eight listed rows excluding "2, additional" are 76+340+47+72+36+29+10+56 = 666, so the statement is in fact correct if "first eight" means the eight main categories). No correction needed there.

## Useful suggestions kept

- A per-unitig table with independent copy-number truth (Hi-C contacts or same-animal ONT) for all 76 category-1 gaps.
- An independent assembler (HiCanu or LJA) at the gap sites to see whether another graph keeps the second copy; Flye collapses haplotypes and is diagnostic only.
- RagTag `patch --fill-only`, run per haplotype, as the established analogue of our read-verified repair.
