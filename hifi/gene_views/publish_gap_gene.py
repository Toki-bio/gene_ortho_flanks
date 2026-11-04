"""Local step: copy TMED4 and EPS8L1 packages into the Tal repo, write IGV sessions and page cards."""
import os, shutil, zipfile, re
V = "C:/work/gene_ortho_flanks/views/gene_views/"; T = "C:/work/Tal/sicista/genes/"
CFG = {
    "TMED4": dict(refs=["tri_hap1", "tri_hap1_fixed", "tri_hap2"], loci={"tri_hap1": "20000-30100", "tri_hap1_fixed": "20000-32400", "tri_hap2": "20000-32400"}),
    "EPS8L1": dict(refs=["tri_hap1", "tri_hap2"], loci={"tri_hap1": "20000-30100", "tri_hap2": "20000-47000"}),
}


def rd(p):
    d = {}; n = None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}


for g, c in CFG.items():
    s = V + g + "/"; d = T + g + "/"
    os.makedirs(d, exist_ok=True)
    aln = rd(s + g + ".aln.fa"); aln.pop("EPS8L1_bet_primary", None)
    with open(d + g + ".aln.fa", "w", newline="\n") as o:
        for k, v in aln.items():
            o.write(f">{k}\n{v}\n")
    bed = [l for l in open(s + g + ".bed") if l.strip() and not l.startswith("EPS8L1_bet_primary")]
    open(d + g + ".aln.bed", "w", newline="\n").writelines(bed)
    names = [f"{g}_{r}" for r in c["refs"]]
    for n in names:
        shutil.copy(s + n + ".sam", d + n + ".sam")
        q = rd(s + "ref.fa")[n]
        open(d + n + ".fa", "w", newline="\n").write(f">{n}\n{q}\n")
        open(d + n + ".bed", "w", newline="\n").writelines(l for l in bed if l.split("\t")[0] == n)
    shutil.copy(s + "ref.fa", d + "ref.fa")
    open(d + "ref.bed", "w", newline="\n").writelines(l for l in bed if l.split("\t")[0] in names)
    for f in ("reads.bam", "reads.bam.bai"):
        shutil.copy(s + f, d + f)
    locus = " ".join(f"{g}_{r}:{c['loci'][r]}" for r in c["refs"])
    open(d + g + ".igv_session.xml", "w", newline="\n").write(f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="ref.fa" hasGeneTrack="false" hasSequenceTrack="true" locus="{locus}" version="8">
    <Resources>
        <Resource path="reads.bam" type="bam"/>
        <Resource path="ref.bed" type="bed"/>
    </Resources>
    <Panel height="600" name="DataPanel" width="1400">
        <Track attributeKey="reads.bam Coverage" autoScale="true" id="reads.bam_coverage" name="HiFi coverage" showReference="false" visible="true"/>
        <Track attributeKey="reads.bam" id="reads.bam" name="HiFi reads ERR16751478, each sequence mapped alone (MAPQ>=20, aligned >=1 kb)" visible="true"/>
    </Panel>
    <Panel height="60" name="FeaturePanel" width="1400">
        <Track attributeKey="ref.bed" id="ref.bed" name="exons, N run, fill" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.8"/>
</Session>
''')
    cwd = os.getcwd(); os.chdir(d)
    with zipfile.ZipFile(g + "_igv.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in [g + ".igv_session.xml", "ref.fa", "ref.bed", "reads.bam", "reads.bam.bai"]:
            z.write(f)
    os.chdir(cwd)

CARDS = '''<section class="card">
  <h2>TMED4 <span class="tag tag-ok">gene on both haplotypes</span></h2>
  <p>The hap1 N run <b>OZ418347.1:96,007,475-96,007,574</b> sits inside the gene and hides 4 of its 5 exons. Hap2 has all 5. The gene-search version fills the N with 2,331 bp copied from hap2 (<b>OZ418087.1:95,802,105-95,804,435</b>), which restores 5 of 5 exons. That fill is the hap2 allele, not hap1's own sequence.</p>
  <p>One read shows what hap1 really has: <b>ERR16751478.985790</b> matches the hap1 flanks (about 11 differences, against 260 on hap2) and crosses the gap. It carries all 5 exons, so the gene is on hap1 too. Between the flanks it has about 2,500 bp, about 170 bp more than the hap2 fill, with 30 mismatches in the aligned part. Hap1 is under-read here: of 13 reads near the gap only 3 are hap1 reads, and the other two end within 10 bp of the gap. One read is below the 3-read repair rule, so hap1's exact sequence is shown but not accepted as a repair.</p>
  <p><b>1. Alignment</b>: trizona hap1, hap2, hap1 fixed, betulina hap1, hap2, and the hap1 read.<br>
  <a class="btn" href="#" onclick="openMSA('sicista/genes/TMED4/TMED4.aln.fa','TMED4 alleles','sicista/genes/TMED4/TMED4.aln.bed')">Open alignment in ViewAlign</a></p>
  <p><b>2. IGV session</b>: before (hap1 with N), after (fixed) and hap2 side by side. Unzip, then File &rarr; Open Session &rarr; TMED4.igv_session.xml.<br>
  <a class="btn green" href="sicista/genes/TMED4/TMED4_igv.zip">Download IGV session (zip)</a></p>
  <p><b>3. HiFi reads in ViewAlign (Reads mode)</b>:<br>
  <a class="btn" href="#" onclick="openRef('TMED4','TMED4_tri_hap1','TMED4 hap1 before, with N')">Before: hap1 with N (53)</a>
  <a class="btn" href="#" onclick="openRef('TMED4','TMED4_tri_hap1_fixed','TMED4 hap1 after, fill from hap2')">After: fixed hap1 (53)</a>
  <a class="btn" href="#" onclick="openRef('TMED4','TMED4_tri_hap2','TMED4 hap2')">Hap2 (60)</a></p>
  <p class="caption">Windows: hap1 OZ418347.1:95,982,475-96,032,574; hap2 OZ418087.1:95,777,189-95,829,000; betulina h1tg000303l:3,583,161-3,637,365 and h2tg001824l:453,931-507,671 (reverse strand). In the fixed view 7 and 8 reads cross the two fill edges, but 4 of the 5 reads crossing the original gap are hap2 reads, so those counts do not prove the fill belongs to hap1.</p>
</section>
<section class="card">
  <h2>EPS8L1 <span class="tag tag-candidate">not repaired, hap1 unknown</span></h2>
  <p>The hap1 N run <b>OZ418348.1:12,671,658-12,671,757</b> stands where hap2 has 16.3 kb (<b>OZ418088.1:13,047,934-13,064,230</b>) holding 18 of the 20 exons. Hap1 shows 2 of 20, hap2 shows 20 of 20. It stays N: the region was classed as tandem and no read crosses it. All 14 reads near the gap fit hap2 better, so I have no read showing what hap1 has there. The gene is present on hap2, and whether hap1 has it is unknown.</p>
  <p><b>1. Alignment</b>: trizona hap1, hap2, betulina hap2. There is no betulina hap1 row: I could not find this region in that assembly beyond short fragments, and the betulina primary assembly is identical to hap2 here (74,398 of 74,398 bp).<br>
  <a class="btn" href="#" onclick="openMSA('sicista/genes/EPS8L1/EPS8L1.aln.fa','EPS8L1 alleles','sicista/genes/EPS8L1/EPS8L1.aln.bed')">Open alignment in ViewAlign</a></p>
  <p><b>2. IGV session</b>: hap1 with N and hap2 side by side. There is no after version. Unzip, then File &rarr; Open Session &rarr; EPS8L1.igv_session.xml.<br>
  <a class="btn green" href="sicista/genes/EPS8L1/EPS8L1_igv.zip">Download IGV session (zip)</a></p>
  <p><b>3. HiFi reads in ViewAlign (Reads mode)</b>:<br>
  <a class="btn" href="#" onclick="openRef('EPS8L1','EPS8L1_tri_hap1','EPS8L1 hap1 with N')">Hap1 with N (46)</a>
  <a class="btn" href="#" onclick="openRef('EPS8L1','EPS8L1_tri_hap2','EPS8L1 hap2')">Hap2 (51)</a></p>
  <p class="caption">Windows: hap1 OZ418348.1:12,646,658-12,696,757; hap2 OZ418088.1:13,022,100-13,093,388; betulina h2tg000165l:1,526,673-1,601,070 (reverse strand). Reads for both genes come from the fished gap reads plus the HiCanu BAM, not a scan of all reads, so a hap1 read could have been missed.</p>
</section>
'''
p = "C:/work/Tal/sicista_genes.html"
t = open(p, encoding="utf-8").read()
if "TMED4" not in t:
    a = t.index('<section class="card">')
    t = t[:a] + CARDS + t[a:]
    t = t.replace("function openReads(hap, title) {", """function openRef(gene, name, title) {
  const b = RAW + 'sicista/genes/' + gene + '/' + name;
  window.open(MSA + '?url=' + encodeURIComponent(b + '.fa') + '&reads=' + encodeURIComponent(b + '.sam') + '&bed=' + encodeURIComponent(b + '.bed') + '&title=' + encodeURIComponent(title), '_blank');
}
function openReads(hap, title) {""")
    open(p, "w", encoding="utf-8", newline="\n").write(t)
for f in ("C:/work/Tal/sicista.html", "C:/work/Tal/sicista/build_page.py"):
    s = open(f, encoding="utf-8").read()
    if "sicista_genes.html" not in s:
        m = re.search(r'<a href="sicista_qc\.html"[^>]*>[^<]*</a>', s)
        print(f, "link anchor", bool(m))
        if m:
            s = s[:m.end()] + ' &middot; <a href="sicista_genes.html" style="color:#fff;">Gene gaps, gene by gene</a>' + s[m.end():]
            open(f, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
