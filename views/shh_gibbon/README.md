# Gibbon SHH 205-bp deletion (Wang ... Wu, Cell 2026, doi 10.1016/j.cell.2025.10.016) - in Sicista?

2026-11-21. The paper's coordinates were not available (full text closed). Location reconstructed:
- UCSC hg38 netNomLeu3 (gibbon Nomascus leucogenys nomLeu3) has a gap of exactly 205 bp at hg38 chr7:155801420-155801625
  (0 bp in gibbon, chr13:106561854). netPanTro6 and netPonAbe3 have no gap there. The only 205-bp gibbon-specific gap in SHH.
- It lies in SHH exon 3 (hg38 155799979-155803726, minus strand), 3' UTR (CDS ends at 155802899).
- netMm10: mouse has the segment (one 10-bp gap inside), mm10 chr5 ~28456480-28456690.
=> this is my reconstruction, not the paper's published coordinates.

Server: therioserver:/home/toki/Sicista2026/shh_gibbon/ (map_shh.sh, step2.sh).
- Query: human hg38 chr7:155799401-155803600 (deletion = offsets 2021-2225) + mouse mm10 chr5:28454401-28458800.
- minimap2 asm20: one MAPQ-60 hit to the SHH coding part in jerboa, Zapus, betulina hap1, betulina primary, trizona hap1/hap2.
  Betulina hap2: no hit (SHH not on a hap2 contig); betulina primary used instead.
- Windows = human offsets 1300-3000 projected through the hit; del.aln.fa = MAFFT L-INS-i of human 1500-2800, mouse 1400-2950, 6 dipodoids.

Result (human columns 2021-2225, the gibbon-deleted segment):
| row | bases aligned in the 205 | identity to human | to mouse |
| mouse | 169 | 53% | - |
| jerboa | 175 | 59% | 55% |
| Zapus | 179 | 60% | 56% |
| Sbet hap1 / primary | 163 | 61% | 60% |
| Stri hap1 / hap2 | 163 | 62% | 60% |
Flanks: right flank 74-79% to human in all dipodoids.
=> the 205 bp are present in both Sicista species (and jerboa, Zapus). The gibbon deletion is gibbon-specific; no matching deletion in Dipodoidea.

ViewAlign files: `shh_3utr.aln.fa` (rows human, mouse, jerboa, Zapus, Sbet hap1, Sbet primary, Stri hap1, Stri hap2; human = hg38
chr7:155800900-155802200 in query offsets 1500-2800) and `shh_3utr.bed` (the 205 bp on the human row, ungapped 521-726).
