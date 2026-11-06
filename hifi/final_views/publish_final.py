"""Local step: copy the views into the Tal repo and write sicista_recovery.html from views.tsv and the final tables.
Style block is read from sicista_genes.html. Every count on the page comes from these tables."""
import csv, os, re, shutil, tarfile, zipfile, collections, html

V = "C:/work/gene_ortho_flanks/views/final_2026-11-06/"
H = "C:/work/gene_ortho_flanks/hifi/final_2026-11-06/"
T = "C:/work/Tal/"
D = T + "sicista/recovery/"


def tsv(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


views = tsv(V + "views.tsv")
gaps = tsv(H + "gap_status.tsv")
fills = {(r["assembly"], r["repair_gap"]): r for r in tsv(H + "new_fills.tsv")}
genes = {r["gene"]: r for r in tsv(H + "final.tsv")}
hapn = {"dep1": "hap1", "dep2": "hap2"}

# data files
os.makedirs(D, exist_ok=True)
for r in views:
    i = r["id"]; s = V + "out/" + i + "/"; d = D + i + "/"
    os.makedirs(d, exist_ok=True)
    for f in ("aln.fa", "view.fa", "view.bam", "other.fa", "other.bam"):
        if os.path.exists(s + f):
            shutil.copy(s + f, d + f)
    for f in ("aln.bed", "view.bed", "other.bed"):
        if os.path.exists(s + f):
            t = open(s + f, encoding="utf-8").read().replace("fill, 0 bp", "joined here, no fill")
            open(d + f, "w", encoding="utf-8", newline="\n").write(t)
    if os.path.exists(s + "igv.tgz"):
        with tarfile.open(s + "igv.tgz") as tf, zipfile.ZipFile(d + i + "_igv.zip", "w", zipfile.ZIP_DEFLATED) as z:
            for m in tf.getmembers():
                if m.isfile():
                    z.writestr(os.path.basename(m.name), tf.extractfile(m).read())


def pos(r):
    return f"{r['contig']}:{int(r['start1']):,}"


def btn(label, js, cls=""):
    return f'<a class="btn{cls}" href="#" onclick="{html.escape(js, quote=True)}">{label}</a>'


def al(r):
    return btn("Alignment", f"aln('{r['id']}')")


def rd(r, f="view", label="Reads"):
    return btn(label, f"reads('{r['id']}','{f}')")


VN = {"own_or_identical": ("own allele or identical", "tag-ok"), "wrong_allele": ("wrong allele", "tag-candidate"),
      "wrong_allele_structure": ("wrong allele", "tag-candidate"), "undecided": ("undecided", "tag-partial"),
      "no_target_cross": ("no own read crosses", "tag-no"), "flanks_tie": ("flanks tie", "tag-no"), "no_phasable_reads": ("no phasable reads", "tag-no")}
ORDER = ["wrong_allele", "own_or_identical", "undecided", "no_target_cross", "flanks_tie", "no_phasable_reads"]
ACT = {"insert": "insert", "join": "join", "genesearch": "gene-search"}
CAUSE = {"tandem": "Tandem array", "dropout": "One parent under-read", "high_copy": "High copy", "fill_repeat": "Fill is a repeat",
         "other": "Other", "ga_tract": "GA tract", "het": "Two alleles unresolved"}
state = {(g["assembly"], g["gap"]): g["final_state"] for g in gaps}

new = [r for r in views if "new" in r["kind"]]
assert len(new) == 78
new.sort(key=lambda r: (ORDER.index(r["verdict_final"]), r["asm"], r["gap"]))
vc = collections.Counter(r["verdict_final"] for r in new)
rows = []
for r in new:
    f = fills[(r["asm"], r["gap"])]
    n, c = VN[r["verdict_final"]]
    strict = "yes" if state[(r["asm"], r["gap"])] == "closed_strict" else "no"
    sub = "" if r["T_x_median"] == "NA" else r["T_x_median"]
    osub = "" if r["O_x_median"] == "NA" else r["O_x_median"]
    rows.append(f'<tr><td>{r["gap"]}</td><td>{hapn[r["asm"]]}</td><td>{pos(r)}</td><td>{ACT[r["action"]]}</td><td class="num">{int(f["fill_len"]):,}</td>'
                f'<td class="num">{f["support"]}</td><td class="num">{r["T"]} / {r["O"]} / {r["tie"]}</td><td class="num">{r["T_cross"]} / {r["O_cross"]}</td>'
                f'<td class="num">{sub}</td><td class="num">{osub}</td><td><span class="tag {c}">{n}</span></td><td>{strict}</td><td>{al(r)} {rd(r)}</td></tr>')
FILLS = "\n".join(rows)
changed = [r for r in new if r["verdict_final"] != r["verdict_fresh"]]
same = 78 - len(changed)
ch = "; ".join(f'{x["gap"]} moved from {VN[x["verdict_final"]][0]} to {VN[x["verdict_fresh"]][0]}' for x in changed)

# open gene gaps
rows = []
gv = [r for r in views if "gene" in r["kind"] and r["state"] == "N"]
for r in sorted(gv, key=lambda x: x["gap"]):
    ex = []
    for g in r["label"].split(","):
        x = genes[g]
        ex.append(f'{g.split("_")[0]} {x["hap1_exons"]} / {x["hap2_exons"]}')
    rows.append(f'<tr><td>{r["gap"]}</td><td>{pos(r)}</td><td>{"<br>".join(ex)}</td><td>{CAUSE[r["lr_cause"]]}</td>'
                f'<td class="num">{r["T"]} / {r["O"]} / {r["tie"]}</td><td class="num">{r["T_cross"]} / {r["O_cross"]}</td>'
                f'<td>{al(r)} {rd(r, "view", "Hap1 reads")} {rd(r, "other", "Hap2 reads")} <a class="btn green" href="sicista/recovery/{r["id"]}/{r["id"]}_igv.zip">IGV</a></td></tr>')
GENES = "\n".join(rows)
ngene = sum(len(r["label"].split(",")) for r in gv)
cross_own = sum(int(r["T_cross"]) for r in gv)
cross_oth = sum(int(r["O_cross"]) for r in gv)

# blocked
cnt = collections.Counter((g["lr_cause"], g["assembly"]) for g in gaps if g["final_state"] == "N")
bl = [r for r in views if r["kind"] == "blocked"]
rows = []
order = sorted({c for c, _ in cnt}, key=lambda c: -(cnt[(c, "dep1")] + cnt[(c, "dep2")]))
for c in order:
    ex = sorted((r for r in bl if r["lr_cause"] == c), key=lambda r: r["asm"])
    for k, r in enumerate(ex):
        head = f'<td rowspan="{len(ex)}">{CAUSE[c]}</td><td class="num" rowspan="{len(ex)}">{cnt[(c, "dep1")]} / {cnt[(c, "dep2")]}</td>' if k == 0 else ""
        rows.append(f'<tr>{head}<td>{r["gap"]}</td><td>{hapn[r["asm"]]}</td><td>{pos(r)}</td><td class="num">{r["T"]} / {r["O"]} / {r["tie"]}</td>'
                    f'<td class="num">{r["T_cross"]} / {r["O_cross"]}</td><td>{al(r)} {rd(r)}</td></tr>')
BLOCK = "\n".join(rows)

nos = next(r for r in views if r["gap"] == "gap0494" and r["asm"] == "dep1")
tm = next(r for r in views if "control" in r["kind"])

style = re.search(r"<style>.*?</style>", open(T + "sicista_genes.html", encoding="utf-8").read(), re.S).group(0)
page = f'''<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sicista trizona final gap recovery</title>
{style}
</head><body>
<header>
  <h1>Sicista trizona final gap recovery</h1>
  <div class="sub">What the last round closed, what it did to genes, and what stays blocked, with alignments and reads for each case. 6 November 2026. Unpublished, for review.
  <a href="sicista_genes.html" style="color:#fff;">&larr; Gene gaps, gene by gene</a></div>
</header>
<main>
<section class="card">
  <h2>Result</h2>
  <table class="tbl"><thead><tr><th>N runs</th><th class="num">Deposited</th><th class="num">Round 2 strict</th><th class="num">Final strict</th><th class="num">Round 2 gene-search</th><th class="num">Final gene-search</th></tr></thead><tbody>
  <tr><td>hap1</td><td class="num">731</td><td class="num">619</td><td class="num"><b>597</b></td><td class="num">334</td><td class="num"><b>292</b></td></tr>
  <tr><td>hap2</td><td class="num">702</td><td class="num">606</td><td class="num"><b>588</b></td><td class="num">319</td><td class="num"><b>286</b></td></tr>
  </tbody></table>
  <p>The last round ran <i>sineclose lr</i> on every gap of both haplotypes, with all assemblies as fill sources and a search of the whole HiFi file. It found 78 fills that round 2 did not have. One gene changed: NOSIP became readable on hap1. {ngene} genes at {len(gv)} gaps are still hidden on hap1, and hap2 shows every one of them.</p>
  <p class="caption">Strict holds joins and inserts only. Gene-search also holds one parent's sequence at gaps where the other parent could not be read, so it is for finding genes. Assemblies: therioserver:/home/toki/gof/final/.</p>
</section>
<section class="card">
  <h2>How to read the views</h2>
  <p><b>Alignment</b> rows: the deposited haplotype with its N run, the final version with the fill (when there is one), the other haplotype, betulina primary, and reads that cross the gap (<code>read_hap1_1</code> and so on: up to 3 of this haplotype, 1 of the other). The annotation track marks the N run, the fill and, for genes, the exons.</p>
  <p><b>Reads</b> opens ViewAlign in Reads mode on the final window, or on the deposited window where the gap is still N. Pink and blue are the two strands, not the haplotypes. Hover over a read to see its name, which starts with its phase: <code>hap1|</code>, <code>hap2|</code>, <code>tie|</code> or <code>nophase|</code>, with <code>_cross</code> added when the read crosses the gap, for example <code>hap1_cross|ERR16751478.240853</code>.</p>
  <p><b>Where to look:</b> views open at column 1. The gap or fill starts 5,000 bp into the window, or 25,000 bp for genes.</p>
  <p class="caption">A read is phased by its differences from each haplotype over the same read bases in the flanks. "Crosses" means at least 400 of the 500 bp on each side of the gap are aligned. Reads were fished again from the whole file with words from both haplotypes, so these counts are new and independent of the run that made the fills. A read made only of repeat cannot be fished. Windows reach 5 kb past the gap, 25 kb for genes. Reads have soft clips and qualities removed to keep the files small.</p>
</section>
<section class="card">
  <h2>NOSIP <span class="tag tag-ok">gene on both haplotypes</span></h2>
  <p>The hap1 N run <b>{nos["contig"]}:{int(nos["start1"]):,}-{int(nos["end1"]):,}</b> hides 6 of the 8 exons. The gene-search fill of {int(nos["fill_len"]):,} bp restores all 8, but it is the hap2 allele. All {nos["O_cross"]} hap2 reads crossing the gap match it with at most 1 substitution.</p>
  <p>One hap1 read shows what hap1 really has. <b>ERR16751478.240853</b> fits the hap1 flanks (7 differences, against 132 on hap2), crosses the gap and carries all 8 exons, so the gene is on hap1 too. Against the fill it has 39 substitutions, 20 indels and one 246 bp insertion. One read is below the 3-read rule, so its sequence is shown here and not used as a repair. This is the same result as TMED4.</p>
  <p>{al(nos)} {rd(nos, "view", "Hap1 final reads ({})".format(nos["view_reads"]))} {rd(nos, "other", "Hap2 reads")} <a class="btn green" href="sicista/recovery/{nos["id"]}/{nos["id"]}_igv.zip">IGV session (zip)</a></p>
  <p class="caption">IGV: unzip, then File &rarr; Open Session &rarr; session.xml. The session has not been opened in IGV yet.</p>
</section>
<section class="card">
  <h2>Genes still hidden on hap1</h2>
  <p>{"No hap1 read crosses any of these %d gaps, and only %d hap2 read crosses one." % (len(gv), cross_oth) if cross_own == 0 and cross_oth == 1 else "CHECK"} So the reads do not show what hap1 holds at these genes. Hap2 has them all; on hap1 they are neither shown nor excluded.</p>
  <table class="tbl"><thead><tr><th>Gap</th><th>hap1 position</th><th>Exons hap1 / hap2</th><th>Cause</th><th class="num">Reads hap1 / hap2 / tie</th><th class="num">Crossing hap1 / hap2</th><th>Views</th></tr></thead><tbody>
{GENES}
  </tbody></table>
  <p class="caption">Exon counts are from the gene call on the strict assembly. At gap0086 most reads tie between the haplotypes, so they cannot be assigned. EPS8L1 also has its older card on the <a href="sicista_genes.html">gene page</a>.</p>
</section>
<section class="card">
  <h2>The 78 new fills</h2>
  <p>{vc["wrong_allele"]} wrong allele, {vc["own_or_identical"]} own allele or identical, {vc["undecided"]} undecided, and {vc["no_target_cross"] + vc["flanks_tie"] + vc["no_phasable_reads"]} that cannot be tested. Wrong allele means reads of this haplotype cross the fill and differ from it: the place is right and the sequence belongs to the other parent. Gene-search takes all 78. Strict takes the joins and inserts except the 3 proven wrong alleles.</p>
  <p>Repeating the phase test on the fresh reads gave the same verdict for {same} of 78; {ch}. Both of those mean the fill cannot be tested.</p>
  <div style="overflow-x:auto"><table class="tbl"><thead><tr><th>Gap</th><th>Hap</th><th>Position</th><th>Kind</th><th class="num">Fill bp</th><th class="num">Junction reads</th><th class="num">Reads own / other / tie</th><th class="num">Crossing own / other</th><th class="num">Subs, own reads</th><th class="num">Subs, other reads</th><th>Verdict</th><th>Strict</th><th>Views</th></tr></thead><tbody>
{FILLS}
  </tbody></table></div>
  <p class="caption">Fill bp below zero: the two contig ends overlap and the join removes that many bases. Junction reads: reads supporting the weaker join in the fill run. Subs: median substitutions inside the fill in reads that span it.</p>
</section>
<section class="card">
  <h2>What stays blocked</h2>
  <p>One example per cause and haplotype, taken as the middle gap of that class in genome order, not chosen by eye. Gaps with a few crossing reads are still below the rule of at least 3 reads at each join.</p>
  <table class="tbl"><thead><tr><th>Cause</th><th class="num">Gaps hap1 / hap2</th><th>Example</th><th>Hap</th><th>Position</th><th class="num">Reads own / other / tie</th><th class="num">Crossing own / other</th><th>Views</th></tr></thead><tbody>
{BLOCK}
  </tbody></table>
</section>
<section class="card">
  <h2>Control</h2>
  <p>TMED4 went through the same steps with a known answer: wrong allele, with read ERR16751478.985790 crossing as a hap1 read and carrying 5 of 5 exons. It came back that way: {tm["T"]} hap1 reads, {tm["O"]} hap2 reads, {tm["T_cross"]} hap1 read and {tm["O_cross"]} hap2 reads crossing, 20 substitutions and 7 indels in that read.</p>
  <p>{al(tm)} {rd(tm, "view", "Reads")} <a class="btn green" href="sicista/recovery/{tm["id"]}/{tm["id"]}_igv.zip">IGV session (zip)</a></p>
</section>
</main>
<script>
const RAW = 'https://raw.githubusercontent.com/Toki-bio/Tal/main/sicista/recovery/';
const MSA = 'https://toki-bio.github.io/MSA-viewer/';
function aln(id) {{
  const b = RAW + id + '/';
  window.open(MSA + '?url=' + encodeURIComponent(b + 'aln.fa') + '&bed=' + encodeURIComponent(b + 'aln.bed') + '&title=' + encodeURIComponent(id), '_blank');
}}
function reads(id, f) {{
  const b = RAW + id + '/' + f;
  window.open(MSA + '?url=' + encodeURIComponent(b + '.fa') + '&reads=' + encodeURIComponent(b + '.bam') + '&bed=' + encodeURIComponent(b + '.bed') + '&title=' + encodeURIComponent(id + ' reads'), '_blank');
}}
</script>
</body></html>
'''
open(T + "sicista_recovery.html", "w", encoding="utf-8", newline="\n").write(page)

# links and the old TMED4 sentence on the gene page
p = T + "sicista_genes.html"
t = open(p, encoding="utf-8").read()
t = t.replace("about 170 bp more than the hap2 fill, with 30 mismatches in the aligned part.", "about 170 bp more than the hap2 fill, with 20 substitutions and 7 indels in the aligned part.")
t = t.replace("Hap1 is under-read here: of 13 reads near the gap only 3 are hap1 reads, and the other two end within 10 bp of the gap.", f"Hap1 is under-read here: fished again from the whole read file, {tm['T']} hap1 and {tm['O']} hap2 reads lie near the gap, and only this one hap1 read crosses it.")
if "sicista_recovery.html" not in t:
    t = t.replace('<a href="sicista.html" style="color:#fff;">&larr; Sicista SINE page</a></div>', '<a href="sicista.html" style="color:#fff;">&larr; Sicista SINE page</a> &middot; <a href="sicista_recovery.html" style="color:#fff;">Final gap recovery</a></div>')
open(p, "w", encoding="utf-8", newline="\n").write(t)
print("fills", dict(vc), "gene gaps", len(gv), "genes", ngene, "cross", cross_own, cross_oth, "blocked rows", len(bl), "unchanged", same, ch)
print("page bytes", len(page))
