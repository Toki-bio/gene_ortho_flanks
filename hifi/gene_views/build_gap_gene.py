#!/usr/bin/env python3
"""Per-gene gap package: deposited hap1 window with the N run, deposited hap2 window, both betulina alleles, the fixed
hap1 if one exists; reads pulled from small sources (repair/fished.fq.gz + indexed HiCanu BAM), mapped to each sequence alone;
exons placed by ssearch36 of the exon sequences against the windows. Writes gene_views/GENE/."""
import subprocess, re, os, collections
G = "/home/toki/gof"; E = "/home/toki/miniforge3/envs"
S = f"{E}/sinederella/bin/samtools"; M = f"{E}/hifiasm/bin/minimap2"; SS = f"{E}/sinederella/bin/ssearch36"
MAFFT = f"{E}/sinederella/bin/mafft"
T1 = "/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna"
T2 = "/home/toki/Sicista2026/demog/genomes/tri_hap2.fna"
B1 = "/home/toki/Sicista2026/assembly/hifiasm_out/sicista_hifiasm.p_ctg.fa"
B2 = f"{G}/betulina/gw/bet_hap2.fa"
HC = f"{G}/asm2/hicanu/tri.contigs.fasta"; HB = f"{G}/asm2/hicanu/tri.contigs.bam"
GENES = {
    "TMED4": dict(id="TMED4_20", c1="OZ418347.1", gs=96007475, ge=96007574, c2="OZ418087.1", w2=(95777189, 95829000),
                  fill=(95802105, 95804435), gsfa=f"{G}/repair/dep1.genesearch.fa"),
    "EPS8L1": dict(id="EPS8L1_237", c1="OZ418348.1", gs=12671658, ge=12671757, c2="OZ418088.1", w2=(13022100, 13093388), fill=None),
}
FLANK = 25000


def sh(c, **k):
    return subprocess.run(c, shell=True, capture_output=True, text=True, **k).stdout


def get(fa, reg):
    return "".join(sh(f"{S} faidx {fa} {reg}").split("\n")[1:]).upper()


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def wfa(p, recs):
    with open(p, "w") as o:
        for n, s in recs:
            o.write(f">{n}\n{s}\n")


def main():
    os.chdir(f"{G}/gene_views")
    Q = {}
    for g, d in GENES.items():
        os.makedirs(g, exist_ok=True)
        s = d["gs"] - FLANK; e = d["ge"] + FLANK
        d["w1"] = (s, e)
        d["h1"] = get(T1, f"{d['c1']}:{s}-{e}")
        d["h2"] = get(T2, f"{d['c2']}:{d['w2'][0]}-{d['w2'][1]}")
        d["gap0"] = d["gs"] - s; d["gap1"] = d["ge"] - s + 1
        assert set(d["h1"][d["gap0"]:d["gap1"]]) == {"N"}, g
        print(g, "hap1", len(d["h1"]), "N in window", d["h1"].count("N"), "hap2", len(d["h2"]), "N", d["h2"].count("N"))
        Q[g] = d["h2"]
    wfa("stage/q2.fa", Q.items())
    # betulina alleles
    for tag, fa in (("bet_hap1", B1), ("bet_hap2", B2)):
        o = sh(f"{M} -x asm20 -c --secondary=no -t 24 {fa} stage/q2.fa 2>/dev/null")
        by = collections.defaultdict(lambda: collections.Counter()); rows = collections.defaultdict(list)
        for l in o.splitlines():
            f = l.split("\t")
            by[f[0]][(f[5], f[4])] += int(f[9]); rows[f[0]].append(f)
        for g, d in GENES.items():
            (c, st), m = by[g].most_common(1)[0]
            r = [x for x in rows[g] if x[5] == c and x[4] == st]
            a = min(int(x[7]) for x in r); b = max(int(x[8]) for x in r)
            q = get(fa, f"{c}:{a+1}-{b}")
            d[tag] = rc(q) if st == "-" else q
            d[tag + "_loc"] = f"{c}:{a+1:,}-{b:,} ({st})"
            print(g, tag, d[tag + "_loc"], "matches", m, "other contigs", by[g].most_common(4)[1:])
    # fixed hap1
    for g, d in GENES.items():
        d["seqs"] = [("tri_hap1", d["h1"]), ("tri_hap2", d["h2"])]
        if d["fill"]:
            fill = get(T2, f"{d['c2']}:{d['fill'][0]}-{d['fill'][1]}")
            d["fx"] = d["h1"][:d["gap0"]] + fill + d["h1"][d["gap1"]:]
            d["fill_len"] = len(fill)
            wfa(f"{g}/fx.fa", [("fx", d["fx"])])
            c = sh(f"{S} faidx {d['gsfa']} {d['c1']} > {g}/gsc.fa; {M} -x asm5 -c --secondary=no {g}/gsc.fa {g}/fx.fa 2>/dev/null; rm -f {g}/gsc.fa")
            print(g, "fixed vs gene-search assembly:", [x.split("\t")[2:4] + x.split("\t")[7:11] for x in c.splitlines()])
            d["seqs"].append(("tri_hap1_fixed", d["fx"]))
    # reads: fished reads + reads the HiCanu contigs used here
    wfa("stage/qh.fa", [(f"{g}_{n}", s.replace("N", "")) for g, d in GENES.items() for n, s in d["seqs"][:2]])
    o = sh(f"{M} -x asm10 -c --secondary=no -t 24 {HC} stage/qh.fa 2>/dev/null")
    reg = sorted({f"{f[5]}:{int(f[7])+1}-{f[8]}" for f in (l.split('\t') for l in o.splitlines()) if int(f_len(f)) >= 2000})
    print("HiCanu regions", reg)
    sh(f"{S} view -b -F 0x900 {HB} {' '.join(reg)} | {S} fastq - 2>/dev/null > stage/hc.fq")
    seen = set(); n = 0
    with open("stage/all.fq", "w") as out:
        for src in (f"zcat {G}/repair/fished.fq.gz", "cat stage/hc.fq"):
            p = subprocess.Popen(src, shell=True, stdout=subprocess.PIPE, text=True)
            while True:
                h = p.stdout.readline()
                if not h:
                    break
                rec = h + p.stdout.readline() + p.stdout.readline() + p.stdout.readline()
                k = h.split()[0]
                if k not in seen:
                    seen.add(k); out.write(rec); n += 1
    print("reads pool", n)
    # exons
    q = {}; nm = None
    for l in open(f"{G}/geneloss2/prep/query.fa"):
        if l[0] == ">":
            nm = l[1:].strip(); q[nm] = []
        else:
            q[nm].append(l.strip())
    ex = dict(l.rstrip("\n").split("\t")[:2] for l in open(f"{G}/geneloss2/prep/exons.tsv"))
    for g, d in GENES.items():
        os.chdir(f"{G}/gene_views/{g}")
        rows = list(d["seqs"]) + [("bet_hap1", d["bet_hap1"]), ("bet_hap2", d["bet_hap2"])]
        names = [(f"{g}_{n}", s) for n, s in rows]
        wfa("rows.fa", names)
        ref = "".join(q[f"{d['id']}|reg|10000"]).upper()
        spans = [tuple(map(int, x.split("-"))) for x in ex[d["id"]].split(",")]
        wfa("exons.fa", [(f"exon{i+1}", ref[a:b]) for i, (a, b) in enumerate(spans)])
        hits = collections.defaultdict(dict)
        for l in sh(f"{SS} -n -m 8 -E 1 -T 8 exons.fa rows.fa 2>/dev/null").splitlines():
            f = l.split("\t")
            if len(f) < 12 or l[0] == "#":
                continue
            e, r = f[0], f[1]; a, b = sorted((int(f[8]), int(f[9]))); bits = float(f[11])
            if e not in hits[r] or bits > hits[r][e][3]:
                hits[r][e] = (a, b, float(f[2]), bits, int(f[7]) - int(f[6]), float(f[10]))
        L = {f"exon{i+1}": b - a for i, (a, b) in enumerate(spans)}
        print(g, "exon lengths", L)
        with open("exons_found.tsv", "w") as o:
            for r, _ in names:
                for e in L:
                    h = hits[r].get(e)
                    ok = h is not None and abs(h[4]) + 1 >= 0.7 * L[e] and h[5] < 1e-3
                    o.write(f"{r}\t{e}\t{L[e]}\t{h[0] if ok else ''}\t{h[1] if ok else ''}\t{h[2] if ok else ''}\n")
                print(" ", r, sum(1 for e in L if (hits[r].get(e) and abs(hits[r][e][4]) + 1 >= 0.7 * L[e] and hits[r][e][5] < 1e-3)), "of", len(L))
        for r, s in names:
            wfa(f"{r}.fa", [(r, s)])
        # reads, each reference alone
        bams = []
        for n, s in d["seqs"]:
            r = f"{g}_{n}"
            sh(f"{M} -ax map-hifi --secondary=no -t 16 {r}.fa ../stage/all.fq 2>/dev/null | {S} view -u -F 4 -q 20 -e 'rlen>=1000' - | {S} sort -o {r}.bam -; {S} index {r}.bam; {S} view -h -F 0x800 {r}.bam > {r}.sam")
            bams.append(f"{r}.bam")
            A = []
            for l in sh(f"{S} view -F 0x800 {r}.bam").splitlines():
                f = l.split("\t"); p = int(f[3]) - 1
                A.append((p, p + sum(int(x) for x, y in re.findall(r"(\d+)([MDN=X])", f[5]))))
            if n == "tri_hap1":
                pts = [("left edge", d["gap0"]), ("right edge", d["gap1"])]
            elif n == "tri_hap1_fixed":
                pts = [("fill left", d["gap0"]), ("fill right", d["gap0"] + d["fill_len"])]
            else:
                pts = []
            print(" ", r, "reads", len(A), [(t, sum(1 for a, b in A if a <= x - 500 and b >= x + 500)) for t, x in pts],
                  "reaching gap from left" if n == "tri_hap1" else "",
                  sum(1 for a, b in A if a < d["gap0"] and b >= d["gap0"] - 50) if n == "tri_hap1" else "")
        sh(f"{S} merge -f -o reads.bam {' '.join(bams)}; {S} index reads.bam")
        wfa("ref.fa", [(f"{g}_{n}", s) for n, s in d["seqs"]])
        sh(f"{MAFFT} --thread 16 --nuc --retree 2 --maxiterate 0 rows.fa > {g}.aln.fa 2> mafft.log")
        with open("info.txt", "w") as o:
            o.write(f"hap1 {d['c1']}:{d['w1'][0]:,}-{d['w1'][1]:,} N {d['gs']:,}-{d['ge']:,} slice {d['gap0']}-{d['gap1']}\n")
            o.write(f"hap2 {d['c2']}:{d['w2'][0]:,}-{d['w2'][1]:,}\nbet_hap1 {d['bet_hap1_loc']}\nbet_hap2 {d['bet_hap2_loc']}\n")
            if d["fill"]:
                o.write(f"fill {d['c2']}:{d['fill'][0]:,}-{d['fill'][1]:,} len {d['fill_len']}\n")
        print(g, "aln bytes", os.path.getsize(f"{g}.aln.fa"))
        os.chdir(f"{G}/gene_views")


def f_len(f):
    return f[10]


main()
print("finished")
