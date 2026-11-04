import sys,collections,re
# input: parsed tsv (name, qs-qe, len, subject#class, id, s-range, bits). Output: deduped calls + class summary
rows=[l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
by=collections.defaultdict(list)
for r in rows:
    qs,qe=map(int,r[1].split("-")); by[r[0]].append((qs,qe,r[3],float(r[4].rstrip("%")),r[5],float(r[6])))
def fam(s):
    n,c=s.split("#"); top=c.split("/")[0]
    if top.startswith("SINE"):
        if re.match(r"(SB_)?(DIP|Dip)",n) or "Dip" in c: return "SINE Dip (tRNA-derived)"
        return {"SINE1":"SINE 7SL (B1/Alu type)","SINE2":"SINE tRNA (B2/ID/B4 type)"}.get(top,top)
    return top
out=[];tot=collections.Counter();cnt=collections.Counter()
for nm,h in by.items():
    h.sort(key=lambda x:-x[5]); acc=[]
    for x in h:
        if any(min(x[1],a[1])-max(x[0],a[0])>0.5*(x[1]-x[0]+1) for a in acc): continue
        acc.append(x)
    for x in sorted(acc):
        out.append((nm,)+x); tot[fam(x[2])]+=x[1]-x[0]+1; cnt[fam(x[2])]+=1
with open(sys.argv[2],"w") as o:
    o.write("indel\tstart\tend\tlen\tbest_match\tclass\t%id\tmatch_range\tbits\n")
    for nm,qs,qe,s,idn,sr,b in out: o.write(f"{nm}\t{qs}\t{qe}\t{qe-qs+1}\t{s.split('#')[0]}\t{fam(s)}\t{idn:.0f}\t{sr}\t{b:.0f}\n")
print(len(out),"elements")
for c,v in tot.most_common(): print(f"{c:22s} n={cnt[c]:3d} bp={v}")
