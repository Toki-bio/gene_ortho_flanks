import csv,collections,re,zipfile,os,shutil
SP=os.getcwd()
V="C:/work/gene_ortho_flanks/views/gene_views/TMEM114/"; D="C:/work/Tal/sicista/genes/TMEM114/"
def fa(p): return ''.join(l.strip() for l in open(p) if l[0]!='>').upper()
h2=fa(V+"hap2.fa"); h1=fa(V+"hap1.fa")
H1="tri_hap1_OZ418350.1_110982846-111075262"; H2="tri_hap2_OZ418090.1_110244645-110400691"
b2=""
for l in open("C:/work/gene_ortho_flanks/views/het_genes/TMEM114.bed"):
    f=l.rstrip("\n").split("\t")
    if f[0]=="TMEM114_tri_hap2" and "exon" in f[3]:
        s,e=int(f[1])+30000,int(f[2])+30000; b2+=f"{H2}\t{s}\t{e}\t{f[3]}\t0\t+\t{s}\t{e}\t0,0,150\n"
        print(f[3],h2[s:s+12],"..",h2[e-6:e])
assert h2[55982:55985]=="ATG"
b2+=f"{H2}\t45708\t109274\tsegment deleted on hap1 (63,566 bp)\t0\t+\t45708\t109274\t200,120,0\n"
b1=f"{H1}\t45707\t45709\tdeletion junction: 63,566 bp of hap2 sequence absent here\t0\t+\t45707\t45709\t200,0,0\n"
print("junction check hap1",h1[45708-20:45708],"| hap2 left",h2[45708-20:45708],"| hap1 right",h1[45708:45728],"| hap2 right",h2[109274:109294])
for d in (V,D):
    open(d+"hap1.bed","w",newline="\n").write(b1); open(d+"hap2.bed","w",newline="\n").write(b2); open(d+"ref.bed","w",newline="\n").write(b1+b2)
    s=open(V+"TMEM114.igv_session.xml").read().replace(":40000-105000",":40000-115000"); open(d+"TMEM114.igv_session.xml","w",newline="\n").write(s)
shutil.copy("TMEM114_indel_TEs.tsv",D+"TMEM114_indel_TEs.tsv")
os.chdir(D)
with zipfile.ZipFile("TMEM114_igv.zip","w",zipfile.ZIP_DEFLATED) as z:
    for f in ["TMEM114.igv_session.xml","ref.fa","ref.bed","reads.bam","reads.bam.bai"]: z.write(f)
os.chdir(SP)
rows=list(csv.DictReader(open("TMEM114_indel_TEs.tsv"),delimiter="\t"))
def grp(r):
    n,c,sp=r["best_match"],r["class"],r["match_species"]
    if re.search("Hydra|Perkinsus|Salmo|Aplysia|Danio",sp): return "Weak match to a non-mammal element"
    if n.startswith(("DIP","Dip")) : return "Dip SINE"
    if c.startswith("SINE1"): return "B1-type SINE"
    if c.startswith("SINE"): return "Other SINE (ID, B2, B4 and similar)"
    if c.startswith("ERV"): return "ERV / LTR"
    return c
S={}
for k in ("hap1_deletion","tri_vs_bet"):
    n=collections.Counter();b=collections.Counter()
    for r in rows:
        if r["comparison"]==k: g=grp(r);n[g]+=1;b[g]+=int(r["hit_len"])
    S[k]=[(g,n[g],b[g]) for g,_ in b.most_common()]; print(k,sum(b.values()),S[k])
def tbl(x): return '<table class="tbl"><thead><tr><th>class</th><th class="num">hits</th><th class="num">bp</th></tr></thead><tbody>'+"".join(f'<tr><td>{g}</td><td class="num">{n}</td><td class="num">{b:,}</td></tr>' for g,n,b in x)+"</tbody></table>"
p="C:/work/Tal/sicista_genes.html"; t=open(p,encoding="utf-8").read()
a=t.index('<section class="card">'); z=t.index("</main>")
R="sicista/genes/TMEM114/"
new=f'''<section class="card">
  <h2>TMEM114 <span class="tag tag-partial">heterozygous deletion</span></h2>
  <p>Deposited hap1 has lost a 63,566 bp segment that holds the whole gene (all 4 exons). The segment is intact on deposited hap2 at <b>OZ418090.1:110,290,353-110,353,918</b>, and both betulina haplotypes carry it, so this is a deletion on trizona hap1, not an insertion on hap2. The hap1 junction is at <b>OZ418350.1:111,028,553</b>. There is no N within 300 kb of the locus on either haplotype, so nothing here can be repaired: the sequence is absent from that allele. 14 HiFi reads cross the hap1 junction by at least 1 kb on both sides, and 8 and 9 cross the two hap2 breakpoints.</p>
  <p class="caption">Sequences in the alignment: hap1 OZ418350.1:111,012,846-111,045,261, hap2 OZ418090.1:110,274,645-110,370,690, betulina h1tg000312l:6,327,125-6,423,689 and h2tg001095l:629,999-722,000 (reverse strand); each matches its assembly exactly. Coordinates in final.tsv are in the repaired assemblies (+254,422 bp on hap1, +709,928 bp on hap2 here), not the deposited ones.</p>
  <p><b>1. Alignment</b>: trizona hap1, trizona hap2, betulina hap1, betulina hap2. It opens at column 1, so scroll to the gene. The last 7 kb (columns after 103,214) are not homologous: the betulina sequences have ended and the two trizona ends are different sequence.<br>
  <a class="btn" href="#" onclick="openMSA('{R}TMEM114.aln.fa','TMEM114 alleles','{R}TMEM114.aln.bed')">Open alignment in ViewAlign</a></p>
  <p><b>2. IGV session</b>: HiFi reads over the hap1 junction and the hap2 segment side by side. Unzip, then File &rarr; Open Session &rarr; TMEM114.igv_session.xml.<br>
  <a class="btn green" href="{R}TMEM114_igv.zip">Download IGV session (zip, 0.4 MB)</a></p>
  <p><b>3. HiFi reads in ViewAlign (Reads mode)</b>: ERR16751478 reads on a 92 kb (hap1) and 156 kb (hap2) window; MAPQ &ge; 20, aligned &ge; 3 kb. Exons, the deleted segment and the junction are in the annotation track.<br>
  <a class="btn" href="#" onclick="openReads('hap1','TMEM114 hap1 (deleted allele) with HiFi reads')">Hap1 reads (24)</a>
  <a class="btn" href="#" onclick="openReads('hap2','TMEM114 hap2 (gene allele) with HiFi reads')">Hap2 reads (53)</a></p>
  <p class="caption">Reads come from the HiCanu read-to-contig BAM (tig00001945 for hap1, tig00006196 for hap2), remapped to the deposited windows. They open at the start of the window, not on the gene.</p>
</section>
<section class="card">
  <h2>Mobile elements in the TMEM114 indels</h2>
  <p>Indels of 50 bp or more were taken from direct pairwise alignments of the four sequences, cut out, and searched with ssearch36 (E &le; 1e-5) against the local RepBase_filtered bank (12,030 sequences), SINEBase and the five Dip consensuses. Full table with positions: <a href="{R}TMEM114_indel_TEs.tsv">TMEM114_indel_TEs.tsv</a>.</p>
  <div class="species-grid">
    <div class="sp-card"><h3>Segment deleted on hap1 (63,566 bp)</h3>{tbl(S["hap1_deletion"])}</div>
    <div class="sp-card"><h3>Trizona vs betulina indels</h3>{tbl(S["tri_vs_bet"])}</div>
  </div>
  <p class="caption">Only 26 % of the deleted segment matches the library. Its elements are old fragments, and none spans the deletion, so no single element explains it. Matches to Hydra, Perkinsus, salmon, Aplysia and zebrafish elements are not credible classifications. Between species the clearest events are full-length DIP0 SINEs (1-222 of 222) in betulina that trizona lacks. Betulina hap1 alone also carries a 5,065 bp L1 (L1_Mur3 ORF2, 74 %) and a DIP0 SINE that betulina hap2 lacks, so these are polymorphic within betulina.</p>
</section>
'''
open(p,"w",encoding="utf-8",newline="\n").write(t[:a]+new+t[z:])
