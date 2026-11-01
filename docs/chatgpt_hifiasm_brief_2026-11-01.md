# Investigation brief: HiFi assembly gaps and possible copy-number misclassification in hifiasm 0.25.0

(Written by ChatGPT, pasted by the owner on 2026-11-01, saved verbatim; it may contain errors. Verification: chatgpt_hifiasm_brief_verification_2026-11-01.md in this folder.)

## Objective

Investigate a suspected mechanism behind missing sequence in one haplotype of a diploid animal genome assembled with **hifiasm 0.25.0**, commit `ec9a8b2`.

The specific problem is:

> Sequence is biologically present on both homologous chromosomes, but appears in only one assembled haplotype. The corresponding region in the other haplotype is represented by an `N` gap after chromosome scaffolding.

Determine whether local HiFi coverage reduction, particularly in GC-rich regions, can cause hifiasm's Hi-C copy-number inference to classify a genuinely two-copy sequence as single-copy, after which Hi-C phasing assigns it to only one haplotype.

Do not assume this hypothesis is correct. Test it against the actual assembly files and the exact source code.

---

## 1. Assembly and source-code context

- Assembler: hifiasm 0.25.0
- Commit: `ec9a8b2`
- Ploidy: diploid
- Logged homozygous coverage threshold: **25**
- Logged primary-contig coverage lower threshold: **21** (85% of 25)
- Suspected problematic local coverage: **13–15 reads**
- Expected haploid coverage under a simple diploid interpretation: approximately **12.5×**

Relevant upstream source files at the exact commit:

- `Purge_Dups.cpp`: https://github.com/chhylp123/hifiasm/blob/ec9a8b2/Purge_Dups.cpp
- `Overlaps.cpp`: https://github.com/chhylp123/hifiasm/blob/ec9a8b2/Overlaps.cpp
- `hic.cpp`: https://github.com/chhylp123/hifiasm/blob/ec9a8b2/hic.cpp
- `rcut.cpp`: https://github.com/chhylp123/hifiasm/blob/ec9a8b2/rcut.cpp
- `Assembly.cpp`: https://github.com/chhylp123/hifiasm/blob/ec9a8b2/Assembly.cpp

Use the exact version above. Do not substitute current `master` code without checking whether the relevant implementation is identical.

## 2. Current mechanistic hypothesis

The proposed chain of events is:

1. A sequence is genuinely present on both homologues.
2. Local HiFi coverage is reduced, potentially due to GC-dependent sequencing bias or another sequence-specific effect.
3. The local observed coverage approaches the expected coverage of a single haplotype rather than the combined coverage of both homologues.
4. hifiasm's Hi-C copy-number inference (`infer_mmhap_copy()`) may classify the sequence as one-copy rather than two-copy.
5. Hi-C phasing assigns that unitig to one haplotype.
6. The other haplotype receives no corresponding copy.
7. A downstream chromosome scaffolder may represent the missing region as an `N` gap.

The numerical observation is suggestive: 13–15× is close to 12.5×, the nominal haploid coverage implied by a homozygous coverage estimate of 25×.

**Important:** This is not yet proof. We must establish what coverage measure `infer_mmhap_copy()` actually uses and whether the observed local coverage values correspond to its input.

## 3. Relevant code already identified

### 3.1 Global coverage and primary-unitig threshold

In `Purge_Dups.cpp`, the homozygous coverage estimate is calculated and stored. The log reports the estimated threshold.

In `Overlaps.cpp`, the primary-unitig logic includes:

- `recover_atg_cov_min = H * 0.85`
- `if_primary_unitig()`
- a per-read exact-overlap coverage calculation
- criteria involving the fraction of reads and overlap bases meeting the threshold

For this assembly, the logged values are:

- H = 25
- 85% threshold = 21.25, reported/used approximately as 21

This logic may affect unitig classification, but do not assume it is the direct cause of the missing haplotype copy. Investigate it separately from the Hi-C copy-number inference.

### 3.2 Hi-C copy-number inference

In `hic.cpp`, the function `gen_mmhap_t()` calls:

`infer_mmhap_copy()`

The relevant code sets:

- `hom_cov` from the global homozygous coverage estimate
- `het_cov = hom_cov / polyploidy`

For a diploid assembly, this gives an expected haploid coverage of approximately 12.5 when the homozygous coverage is 25.

The code then marks a unitig inferred to have the full ploidy copy number as fully assigned.

The critical function to inspect is:

**`infer_mmhap_copy()` in `Overlaps.cpp`.**

Its complete implementation has not yet been verified. Determine exactly:

1. What data it uses to estimate copy number.
2. Whether it uses raw read depth, overlap coverage, read multiplicity, or another quantity.
3. How it uses the `het_cov` argument.
4. Its thresholds and rounding rules.
5. How it distinguishes one-copy from two-copy unitigs.
6. Whether local coverage variation can change the inferred copy number.
7. Whether the result is affected by unitig length, coverage variance, repeat content, or read overlap quality.

Do not infer its behaviour from the function name or caller alone. Read the complete function body.

### 3.3 Hi-C phasing

In `hic.cpp`, inspect these functions:

- `gen_mmhap_t()`
- `infer_mmhap_copy()`
- `gen_mmhap_bub()`
- `purge_phase()`
- `purge_phase_0()`
- `hic_short_align_mmhap()`

The already identified logic in `purge_phase()` is important.

For the first haplotype, the procedure performs graph partitioning using Hi-C links. For the second haplotype, an early-return branch assigns remaining unassigned unitigs to the second haplotype.

This suggests that a unitig classified as one-copy may be assigned to one haplotype without being duplicated into the other.

However, verify this interpretation against the complete code and data structures. In particular, determine how the inferred copy number (`m`) and assignment count (`n`) control the number of haplotype assignments.

### 3.4 Confirm the actual execution path

It is essential to establish whether the user's assembly actually invokes `hic_short_align_mmhap()`.

Trace the call chain from the relevant command-line options and top-level assembly function to the Hi-C phasing implementation.

Do not assume that every Hi-C assembly mode uses this route. Check whether another or older Hi-C path is selected under any relevant conditions.

---

## 4. Assembly gap classification

The user has classified 731 gaps into the following categories:

| Category | Description | Count |
|---|---|---:|
| 1 | Same sequence on both copies, but represented only once | 76 |
| 2 | One copy cannot be read at a site; coverage fades towards the gap | 340 |
| 2, additional | Similar cases where no assembly crosses the site | 65 |
| 3 | No sequence missing; scaffolder inserted a placeholder | 47 |
| 4 | Tandem repeat array too long or identical to resolve | 72 |
| 5 | Multiple repeat copies collapsed into one | 36 |
| 6 | Male X chromosome / hemizygous sequence | 29 |
| 7 | Long pure GA/TC tract adjacent to gap | 10 |
| 8 | Repeat sequence present in one haplotype but unresolved or absent in the other | 56 |

The first eight category counts sum to 666; category 2's additional 65 locations bring the total to 731. Check whether categories overlap and preserve the user's counting convention.

Category 1 is the primary focus of this investigation.

Category 2 is a related but distinct problem: coverage appears to fade towards a gap, potentially reflecting sequence-specific read dropout or difficulty sequencing the region.

Do not conflate:
- failure to sequence one homologous copy;
- failure to assemble a sequence that was sequenced;
- copy-number misclassification during phasing;
- repeat collapse;
- gaps introduced by downstream scaffolding.

---

## 5. Relevant observations from the assembly

### Category 1

- 76 cases identified.
- 65 were repaired by the user's subsequent procedure.
- Some candidate regions show approximately 13–15 reads where approximately 25× might be expected for two-copy sequence.
- The hypothesis is that locally reduced coverage makes shared sequence appear single-copy.

### Category 2

- Coverage decreases towards the gap.
- An example has estimated k-mer copy number around 0.52 at the gap and 0.80 approximately 9 kb away.
- Candidate causes include expanded STRs, G-quadruplex-associated regions, GA-rich sequence, and other sequence-dependent sequencing difficulties.

### GC observations

The user reports the following local coverage relationships:

- At 36% GC, two-copy DNA has approximately 1.23× average k-mer count.
- At 52% GC, it has approximately 0.73× average k-mer count.

These observations should be checked against the underlying data and the exact definition of the k-mer metric. Do not automatically equate k-mer count with the coverage measure used internally by hifiasm.

---

## 6. Tasks to perform on the server

Start by inspecting the actual files. Do not assume paths, filenames, or file naming conventions.

### Task A — Inventory

Identify:

- hifiasm log files;
- primary and haplotype GFA files;
- Hi-C-related GFA outputs;
- any intermediate unitig or bubble files retained;
- the relevant FASTA files;
- read files or read-alignment files, if available;
- downstream scaffolded assemblies containing the `N` gaps.

Record the exact paths and file sizes. Do not modify or overwrite anything.

### Task B — Identify representative category-1 cases

Find a few representative category-1 regions, preferably including:

- one of the 65 repaired cases;
- one unrepaired case;
- one case with particularly clear 13–15× local coverage.

For each case, establish:

- genomic/scaffold coordinates;
- sequence and orientation;
- which haplotype contains the sequence;
- whether the other haplotype has a corresponding unitig or graph path;
- local GC content;
- local read depth;
- k-mer copy-number evidence;
- read support across the region and its flanks.

Use exact coordinates and sequence identifiers from the files.

### Task C — Trace the graph and unitigs

For each selected case:

1. Identify the relevant unitig(s) in the hifiasm graph.
2. Determine whether the sequence is represented once or twice in the graph.
3. Determine whether the unitig was considered one-copy or two-copy by the Hi-C code, if this information can be recovered.
4. Trace its haplotype assignment.
5. Establish whether the missing haplotype copy was absent before scaffolding or lost during a later step.

Do not infer graph copy number solely from the final FASTA representation.

### Task D — Test the coverage hypothesis

For representative cases, compare local coverage with suitable controls.

Where possible, compare:

- the problematic region;
- nearby unique sequence;
- sequence of similar GC content;
- regions known to be present in both haplotypes;
- regions known to be genuinely single-copy.

Determine whether local coverage reduction is sufficient to explain a one-copy classification under the actual `infer_mmhap_copy()` implementation.

### Task E — Inspect the downstream scaffolder

The final `N` gaps may have been introduced after hifiasm.

Identify the scaffolder and its input files. Establish whether:

- both haplotype sequences were already missing/unequal before scaffolding;
- the scaffolder inserted a fixed-size placeholder;
- a real sequence was present in an earlier assembly stage but omitted later.

Do not attribute scaffolder-created Ns to hifiasm without evidence.

---

## 7. Required source-code verification

Retrieve the exact hifiasm 0.25.0 source and inspect the complete bodies of the relevant functions.

At minimum, provide:

1. Full implementation of `infer_mmhap_copy()`.
2. Its declaration and all call sites.
3. The complete relevant parts of `gen_mmhap_t()`.
4. The complete relevant parts of `purge_phase()` and `purge_phase_0()`.
5. The caller chain leading to `hic_short_align_mmhap()`.
6. The code that writes the relevant haplotype GFA/FASTA outputs.

Use line-numbered source excerpts and cite the exact commit. If the source cannot be retrieved, say so explicitly rather than reconstructing the implementation from partial snippets.

---

## 8. Questions that the investigation must answer

1. Does `infer_mmhap_copy()` infer copy number from a local coverage measure that can be depressed by GC bias?
2. If the true two-copy sequence has local coverage of 13–15× while the global homozygous estimate is 25×, can it be classified as one-copy?
3. What exact threshold or calculation determines this?
4. If a unitig is classified as one-copy, can the Hi-C phasing algorithm assign it to only one haplotype?
5. Can Hi-C contacts rescue a two-copy classification, or do they only determine which haplotype receives a unitig already classified as one-copy?
6. Does the user's actual command-line mode invoke the mmhap route?
7. At what assembly stage does the second copy disappear?
8. Is the `N` gap introduced by hifiasm or by downstream scaffolding?
9. Can the proposed mechanism explain the category-1 cases, and what evidence would falsify it?
10. Are there alternative explanations that fit the observations better?

---

## 9. Evidence and reporting standards

Separate conclusions into three levels:

**Confirmed**
- Directly supported by the exact source code or assembly files.

**Strongly supported**
- Consistent with measured data and source logic, but one link in the causal chain remains unverified.

**Speculative**
- Plausible explanation requiring additional evidence.

For every important conclusion, provide:
- source file and line numbers;
- exact command or analysis used;
- relevant file path and sequence identifiers;
- numerical evidence;
- limitations.

Do not invent paths, file contents, thresholds, function behaviour, or results. Do not change the assembly files.

## Expected deliverable

Produce a concise technical report containing:

1. A verified source-code explanation.
2. A reconstruction of at least two representative category-1 cases.
3. A diagram of the proposed causal chain, with confirmed and unconfirmed steps distinguished.
4. A conclusion on whether local coverage loss can cause erroneous single-copy classification during Hi-C phasing.
5. Alternative explanations and tests that distinguish them.
6. A minimal list of additional analyses needed, if the available files are insufficient.

The goal is to determine what actually happened in this assembly, not merely to find a plausible explanation in the source code.
