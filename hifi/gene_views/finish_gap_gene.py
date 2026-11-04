#!/usr/bin/env python3
"""After build_gap_gene.py: TMED4 hap1-phased read as its own row, EPS8L1 betulina primary row, exon BEDs, MAFFT."""
import subprocess, re, os, collections
G = "/home/toki/gof/gene_views"; E = "/home/toki/miniforge3/envs"
S = f"{E}/sinederella/bin/samtools"; M = f"{E}/hifiasm/bin/minimap2"; SS = f"{E}/sinederella/bin/ssearch36"; MAFFT = f"{E}/sinederella/bin/mafft"
P = "/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa"


def sh(c):
    return subprocess.run(c, shell=True, capture_output=True, text=True).stdout


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def rd(p):
    d = {}; n = None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}


def wfa(p, d):
    with open(p, "w") as o:
        for k, v in d.items():
            o.write(f">{k}\n{v}\n")


# TMED4: the hap1 read
os.chdir(f"{G}/TMED4")
rows = rd("rows.fa")
for l in open("TMED4_tri_hap1.sam"):
    f = l.split("\t")
    if f[0].endswith(".985790"):
        break
seq = f[9]; q = 0; ins = None; pos = int(f[3]) - 1
for a, b in re.findall(r"(\d+)([MIDSH=X])", f[5]):
    a = int(a)
    if b == "I" and a >= 1000:
        ins = (q, q + a, pos)
    if b in "MIS=X":
        q += a
    if b in "MD=X":
        pos += a
print("read", f[0], "flag", f[1], "len", len(seq), "insert in read", ins)
wfa("insert.fa", {"hap1_read_insert": seq[ins[0]:ins[1]]})
rows["TMED4_hap1_read"] = seq
wfa("fill.fa", {"fill": rows["TMED4_tri_hap1_fixed"][25000:25000 + 2331]})
print(sh(f"{M} -x asm20 -c --cs=short fill.fa insert.fa 2>/dev/null | cut -f1-12"))
for l in sh(f"{M} -x asm20 -c --cs=short fill.fa insert.fa 2>/dev/null").splitlines():
    cs = [x for x in l.split("\t") if x.startswith("cs:Z:")][0]
    print("  indels >=20 between fill and hap1 read:", [(o, len(a)) for o, a in re.findall(r"([+\-])([a-z]+)", cs) if len(a) >= 20], "mismatches", cs.count("*"))
wfa("rows.fa", rows)

# EPS8L1: betulina primary replaces the 842 bp scrap
os.chdir(f"{G}/EPS8L1")
rows = rd("rows.fa")
del rows["EPS8L1_bet_hap1"]
x = sh(f"{S} faidx {P} ptg000065l:23595135-23669532")
rows["EPS8L1_bet_primary"] = rc("".join(x.split("\n")[1:]).upper())
wfa("rows.fa", rows)
for f in ("EPS8L1_bet_hap1.fa", "EPS8L1_bet_hap1.bam"):
    if os.path.exists(f):
        os.remove(f)
print("primary vs bet_hap2:", sh(f"awk '/^>EPS8L1_bet_primary/{{p=1;print;next}}/^>/{{p=0}}p' rows.fa > a.fa; awk '/^>EPS8L1_bet_hap2/{{p=1;print;next}}/^>/{{p=0}}p' rows.fa > b.fa; {M} -x asm20 -c --secondary=no b.fa a.fa 2>/dev/null | cut -f3,4,8,9,10,11; rm -f a.fa b.fa"))

# exons, BEDs, alignments
for g in ("TMED4", "EPS8L1"):
    os.chdir(f"{G}/{g}")
    L = {k: len(v) for k, v in rd("exons.fa").items()}
    hits = collections.defaultdict(dict)
    for l in sh(f"{SS} -n -m 8 -E 1 -T 8 exons.fa rows.fa 2>/dev/null").splitlines():
        f = l.split("\t")
        if len(f) < 12 or l[0] == "#":
            continue
        a, b = sorted((int(f[8]), int(f[9]))); bits = float(f[11])
        if f[0] not in hits[f[1]] or bits > hits[f[1]][f[0]][3]:
            hits[f[1]][f[0]] = (a, b, float(f[2]), bits, abs(int(f[7]) - int(f[6])) + 1, float(f[10]))
    bed = []
    rows = rd("rows.fa")
    for r in rows:
        found = []
        for e, n in L.items():
            h = hits[r].get(e)
            if h and h[4] >= 0.7 * n and h[5] < 1e-3:
                found.append(f"{e}({h[2]:.0f}%)")
                bed.append(f"{r}\t{h[0]-1}\t{h[1]}\t{g} {e} ({n} bp)\t0\t+\t{h[0]-1}\t{h[1]}\t0,0,150")
        print(g, r, len(rows[r]), f"{len(found)}/{len(L)}", " ".join(found))
    bed.append(f"{g}_tri_hap1\t25000\t25100\tN run, 100 bp\t0\t+\t25000\t25100\t200,0,0")
    if g == "TMED4":
        bed.append("TMED4_tri_hap1_fixed\t25000\t27331\tfill copied from hap2, 2,331 bp\t0\t+\t25000\t27331\t200,120,0")
    open(f"{g}.bed", "w").write("\n".join(bed) + "\n")
    sh(f"{MAFFT} --thread 16 --nuc --retree 2 --maxiterate 0 rows.fa > {g}.aln.fa 2> mafft.log")
    print(g, "aln", {k: len(v) for k, v in rd(f"{g}.aln.fa").items()})
    print(sh(f"{S} idxstats reads.bam"))
print("finished")
