#!/usr/bin/env python3
"""Windows and fishing words for sv_window.py on the asm25 haplotypes.
Region windows: every gap-fill (FILLS_vs_HAP.paf, query aligned end to end, MAPQ >= 20) +- PAD on that haplotype,
merged when they overlap. Controls: N_CTRL random windows of 2*PAD on contigs > 5 Mb of hap1, away from the regions.
Fishing words: 25-mers every STEP bp of each window, counted with jellyfish --if in the window's own haplotype; kept
when unique there. Writes OUT/windows.bed, OUT/windows.fa, OUT/anchors.txt.
Usage: sv_prep.py OUT_DIR hap1=FA:PAF hap2=FA:PAF"""
import os
import random
import subprocess
import sys

PAD, N_CTRL, STEP, K = 300000, 4, 50, 25
out = sys.argv[1]
haps = {a.split("=")[0]: a.split("=")[1].split(":") for a in sys.argv[2:]}
os.makedirs(out, exist_ok=True)
comp = str.maketrans("ACGT", "TGCA")


def fai(fa):
    return {l.split()[0]: int(l.split()[1]) for l in open(fa + ".fai")}


wins = []
for h, (fa, paf) in haps.items():
    ln = fai(fa)
    for l in open(paf):
        f = l.split("\t")
        if f[0] == "mid" or int(f[11]) < 20 or int(f[2]) > 100 or int(f[3]) < int(f[1]) - 100:
            continue
        t, ts, te = f[5], int(f[7]), int(f[8])
        wins.append([h, t, max(0, ts - PAD), min(ln[t], te + PAD), f[0].split("__")[0]])
wins.sort()
merged = []
for w in wins:
    if merged and merged[-1][0] == w[0] and merged[-1][1] == w[1] and w[2] <= merged[-1][3]:
        merged[-1][3] = max(merged[-1][3], w[3]); merged[-1][4] += "+" + w[4]
    else:
        merged.append(w)
random.seed(3)
ln = fai(haps["hap1"][0])
big = [c for c, n in ln.items() if n > 5_000_000]
while sum(1 for w in merged if w[4].startswith("CTRL")) < N_CTRL:
    c = random.choice(big)
    s = random.randrange(0, ln[c] - 2 * PAD)
    if not any(w[1] == c and w[2] < s + 2 * PAD and s < w[3] for w in merged):
        merged.append(["hap1", c, s, s + 2 * PAD, f"CTRL{sum(1 for w in merged if w[4].startswith('CTRL')) + 1}"])
with open(os.path.join(out, "windows.bed"), "w") as fh:
    for h, c, s, e, lab in merged:
        fh.write(f"{h}\t{c}\t{s}\t{e}\t{lab}\n")
anchors = set()
with open(os.path.join(out, "windows.fa"), "w") as wf:
    for h in haps:
        ws = [w for w in merged if w[0] == h]
        if not ws:
            continue
        reg = [f"{c}:{s + 1}-{e}" for _, c, s, e, _ in ws]
        seqs = subprocess.run(["samtools", "faidx", haps[h][0]] + reg, capture_output=True, text=True, check=True).stdout
        cand, name, buf = set(), None, []
        recs = []
        for l in seqs.splitlines() + [">"]:
            if l.startswith(">"):
                if name:
                    recs.append((name, "".join(buf).upper()))
                name, buf = l[1:], []
            else:
                buf.append(l)
        for (name, s), w in zip(recs, ws):
            wf.write(f">{w[4]}|{h}|{name}\n{s}\n")
            for i in range(0, len(s) - K, STEP):
                k = s[i:i + K]
                if "N" not in k and len(set(k)) > 2:
                    cand.add(min(k, k.translate(comp)[::-1]))
        cf = os.path.join(out, f"cand_{h}.fa")
        with open(cf, "w") as fh:
            fh.write("".join(f">{i}\n{k}\n" for i, k in enumerate(sorted(cand))))
        jf = os.path.join(out, f"cand_{h}.jf")
        subprocess.run(["jellyfish", "count", "-m", str(K), "-C", "-s", str(4 * len(cand)), "-t", "16", "--if", cf,
                        "-o", jf, haps[h][0]], check=True)
        dump = subprocess.run(["jellyfish", "dump", "-c", jf], capture_output=True, text=True, check=True).stdout
        anchors |= {k for k, n in (x.split() for x in dump.splitlines()) if n == "1"}
        os.remove(jf)
        print(f"{h}: {len(ws)} windows, {len(cand)} candidate words, unique so far {len(anchors)}")
with open(os.path.join(out, "anchors.txt"), "w") as fh:
    fh.write("\n".join(sorted(anchors)) + "\n")
for w in merged:
    print("\t".join(map(str, w)))
