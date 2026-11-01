# hifiasm_gc: GC correction of the "primary unitig" coverage test inside hifiasm 0.25.0

Patch `hifiasm_gc_corr.diff` against hifiasm 0.25.0 (commit ec9a8b2, `Overlaps.cpp` only).

## What it changes and why

hifiasm keeps a unitig in both haplotypes only if `if_primary_unitig` (Overlaps.cpp ~18895) calls it
"primary": per read, C = sum of exact-overlap lengths / read length; a read qualifies when
C >= recover_atg_cov_min (0.85 x hom coverage = 21 in the S. trizona run). A unitig that is not primary
is deleted from the haplotype whose Hi-C/trio labels it lacks (Overlaps.cpp ~20228 and ~20315).
In this HiFi data two-copy coverage falls with GC (relative to the genome two-copy depth: 36 % GC 1.18-1.23,
44 % 0.94, 52 % 0.73-0.78, 56 % 0.72), so GC-rich homozygous stretches (rd 13-15) fail the test, are placed in
one haplotype only and leave a hole in the other.

The patch divides each read's C by f(GC of that read) before the threshold comparison:

    C' = (C_bases / R_bases) / f(GC)

f is a linear interpolation of a GC table, clamped at the ends. Default table (HIFIASM_GC_CORR=on):
36:1.205, 44:0.94, 52:0.755, 56:0.72. Read GC is computed once (pthread_once, kt_for over all reads) from the
corrected reads held in R_INF, on the first call of `if_primary_unitig`; a log line prints the table and the
distribution of factors.

Switch: environment variable `HIFIASM_GC_CORR`
- unset, empty, `0`, `off`: stock behaviour (integer C, identical code path; log line "off").
- `on` or `1`: default table.
- `"gc:f,gc:f,..."`: own table, GC in %, ascending.

Scope: only `if_primary_unitig` (the gate that decides deletion from a haplotype). `set_r_het_flag`
(Overlaps.cpp ~14016, reads with coverage <= 0.75 H flagged heterozygous, used for bubble identification in
Hi-C phasing and purge_dups h_status) is left unchanged on purpose: (1) the measured placement switch
(rd 16-17) matches the 21 gate, not the het flag; (2) the het flag only says which reads may carry phasing
signal; a GC-rich homozygous unitig wrongly flagged het still gets a Hi-C label, but the primary gate then
protects it from deletion; (3) correcting the het flag too would change the phasing input genome-wide and make
the effect of the one change impossible to attribute. `get_utg_attributes` (commented-out caller) and
`recover_utg_by_coverage` (not on the diploid Hi-C path) are not touched.

## Build (therioserver)

System and conda headers clash with hifiasm's `#define __ASSEMBLY__` (Assembly.h) via
linux/sched/types.h; build with `CPPFLAGS="-include pthread.h"` (see `build.sh`). Stock and patched binaries
are built the same way: therioserver:/home/toki/gof/hifiasm_gc/{stock,gc}/hifiasm.

## Files

- `hifiasm_gc_corr.diff` the patch (v2 = v1 + `HIFIASM_GC_DEBUG`, a read-name list; every
  `if_primary_unitig` call on a unitig holding a listed read prints a `[GCDBG]` line). `hifiasm_gc_corr_v1_no_debug.diff`
  = the binary used for the main run. `build.sh` (stock + v1), `build2.sh` (v2 + toy checks).
- `toy.sh` toy (fished reads of the 10 regions + 4 controls, ends/sv/fished.fa, 9702 reads; first 1 M Hi-C pairs).
- `cpbins.sh`, `run_gc.sh` full run reusing asm25's ec/ovlp/Hi-C bins (copied); `run_dbg.sh` debug rerun.
- `eval.sh`, `pre_asm25.sh`, `winpaf.py`, `gapcmp.py`, `utgcmp.sh`, `readdiff.py` evaluation vs asm25.

Server: therioserver:/home/toki/gof/hifiasm_gc/ (asm_gc = GC run, asm_dbg = debug rerun, eval = evaluation).

## Results (2026-11-01)

Toy: stock vs patched with correction off, all 30 output gfa/bed files byte-identical (v1 and v2); v2 on = v1 on.

Full run (HIFIASM_GC_CORR=on, -t 48, same bins and command as asm25): 4,068,421 reads, mean factor 1.012
(the table is centred on the genome). p_utg identical to asm25; Hi-C partition identical (same mc_solve edge counts);
dedup_base identical (1,402,293). The debug rerun reproduced the run byte for byte.

- Haplotypes: hap1 2,629,762,683 -> 2,629,768,663 bp (+6 kb), hap2 2,782,190,923 -> 2,782,572,016 bp (+381 kb;
  806 -> 804 contigs). Reads moved: hap1 +12/-17, hap2 +129/-19, i.e. two unitigs flipped (utg009345l rd17,
  utg002509l rd18).
- 10 studied regions (windows of the piece-carrying asm25 haplotype, minimap2 asm5): best single-contig coverage in
  the other haplotype unchanged at every region (R01 0.337, R04 0.634, R07 0.426, R00 0.480, R02 0.511, R05 0.653,
  R06 0.412, R09 0.521; R08 and 4 controls 1.0 in both, as before).
- Deposited-gap test (gapfill.py, flanks of every N run): dep1 731 gaps, dep2 702 gaps, status of every gap identical
  between asm25 and the GC run; category 1 (repair action "insert"): dep1 75 = 51 one haplotype / 24 both,
  dep2 65 = 42 / 22 / 1 none, unchanged.
- Unitig rule table (utg_rule_test.py): share placed in both haplotypes identical in every rd x GC bin except rd 17
  (0.74 -> 0.75) and rd 18 GC >= 48 (0.80 -> 0.92). Heterozygous control: all bins rd <= 16 unchanged
  (one-copy unitigs rd 9-13: 0.00-0.04 of length in both haplotypes, as in asm25).
- compleasm (contigs, odb10) and Merqury (both haplotypes, sqc/hifi.meryl) were still running when this was written
  (server load ~200): eval/compleasm.txt, eval/mq_gc/, eval/mq_asm25/. With 0.39 Mb changed out of 5.41 Gb they
  cannot differ by more than that.

Why it does not work (debug run, reads of the R01/R04/R07 piece unitigs, all removed, ret=0):

| piece | reads | mean C (hifiasm) | mean factor | mean C' | reads C >= 21 | reads C' >= 21 | needed (80 %) |
|---|---|---|---|---|---|---|---|
| R01 | 85 | 13.30 | 0.784 | 16.99 | 0 | 15 | 68 |
| R04 | 27 | 14.23 | 0.833 | 17.08 | 0 | 7 | 22 |
| R07 | 38 | 16.89 | 0.793 | 21.30 | 0 | 21 | 31 |

hifiasm's own per-read coverage of these two-copy pieces is 0.51-0.65 of H (25): the measured GC curve explains a
factor ~0.8, the rest (another ~0.7-0.8) is not GC-predictable (allele-specific/local dropout plus hifiasm counting
only exact overlaps). Reaching the 80 % rule would need a gate near 0.5 H, which is where one-copy (heterozygous)
unitigs sit (rd 11-14), so any stronger correction or lower --hom-cov would duplicate heterozygous sequence.

Verdict: in-assembler GC correction of the primary-unitig gate is safe (no measurable false duplication: at most
+0.39 Mb in both haplotypes) but does not fill any of the category-1 holes. Fix them post hoc (read-verified copy of
the fill, repair_gaps.py).
