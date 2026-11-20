# Jerboa Rsrc1/Shox2 deletion: right-end junction in Dipodoidea

Question (2026-11-20): is the ~1 kb deletion in the *Rsrc1* intron ~285 kb upstream of *Shox2*, reported in jerboas and
*Zapus* by Guerrero ... Hiller, Cooper (bioRxiv 2026, doi 10.64898/2026.09.17.752424), also in *Sicista*?

`rightend.aln.fa`: 600 bp around the right end, MAFFT L-INS-i (`--localpair --maxiterate 1000`). Rows, all oriented like mouse:

| row | genome | window |
|---|---|---|
| mouse_mm10 | mm10 | chr3:67267331-67267931 |
| jerboa_Jjac | mJacJac1.mat.Y.cur | NC_059112.1:66958168-66958768 (+) |
| Zapus_hud | mZapHud2 | CM176106.1, window offsets 7356-7956 of 122166488-122182487 (- strand) |
| Sbet_hap1 / Sbet_hap2 | *S. betulina* hifiasm hap1 / hap2 | h1tg001068l / h2tg000110l, window offsets 7262-7862 |
| Stri_hap1 / Stri_hap2 | *S. trizona* mSicTri.1 / deposited hap2 | OZ418350.1 / OZ418090.1, window offsets 7277-7877 / 7267-7867 |

`rightend.bed`: per-row ungapped 0-based marks: the 23-bp block kept in every genome (mouse 6147-6169 of the query
window = mm10 chr3:67267635-67267657); in mouse, the right part of the paper's 1228-bp segment and the ~139 bp
retained piece.

Findings: identity mouse vs every dipodoid is 73-86% from mouse ~6090 rightwards and background (28-48%) leftwards,
while jerboa-*Sicista* stay 67-89%. The right end of the deletion is at mouse window offset ~6090 (mm10 chr3:~67267578),
the same in jerboa, *Zapus* and all four *Sicista* haplotypes. 6228 - 6090 = 138, the paper's retained 139 bp.
*Zapus* has an extra ~100 bp GC-rich stretch 50-170 bp left of the end; MAFFT spreads it across the row.

Full work: `C:\work\hylomys_ont\shox2_deletion\` (README.md, scripts, dot plots); server
`therioserver:/home/toki/Sicista2026/shox2_deletion/`.
