import sys, statistics, re
D='/home/toki/gof/views_het/'
def readfa(p):
    n=[];s=[]
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): n.append(l[1:]);s.append([])
        else: s[-1].append(l)
    return n,[''.join(x).upper() for x in s]
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
cfg={'TMEM114':('/home/toki/gof/hifi/alns/R18.aln.fa',22178,[(0,220),(523,81),(24751,138),(27871,233)]),
     'SMCO1':('/home/toki/gof/hifi/alns/R14.aln.fa',9737,[(0,50),(4039,150),(4843,445)])}
colors={'TMEM114':'0,0,150','SMCO1':'150,0,0'}
for g,(ref,st,ex) in cfg.items():
    rn,rs=readfa(ref); r0=rs[0].replace('-','')
    exons=[r0[st+o:st+o+l] for o,l in ex]
    names,seqs=readfa(D+g+'.aln.raw.fa')
    ung=[s.replace('-','') for s in seqs]
    bed=[]; rep=[]
    for nm,u in zip(names,ung):
        for i,e in enumerate(exons,1):
            for strand,q in (('+',e),('-',rc(e))):
                k=10; votes={}
                for p in range(0,len(q)-k+1):
                    sub=q[p:p+k]; x=u.find(sub)
                    while x!=-1 and len(votes)<5000:
                        votes[x-p]=votes.get(x-p,0)+1; x=u.find(sub,x+1)
                if votes:
                    d,v=max(votes.items(),key=lambda t:t[1])
                    # allow small drift: count matches
                    win=u[d:d+len(q)] if d>=0 else ''
                    m=sum(a==b for a,b in zip(win,q))
                    if len(win)==len(q) and m/len(q)>=0.8 and v>=3:
                        bed.append((nm,d,d+len(q),f'exon{i}',strand,m/len(q)))
                        rep.append(f'{nm} exon{i}({len(q)}bp) {d}-{d+len(q)} {strand} id={m/len(q):.3f}')
                        break
    # gaps
    gaprep=[]
    for nm,s in zip(names,seqs):
        pos=0;
        for m in re.finditer(r'-{300,}',s):
            before=len(s[:m.start()].replace('-',''))
            gaprep.append(f'{nm} gap {m.end()-m.start()} cols at ungapped pos {before} (aln cols {m.start()}-{m.end()})')
    print(g,'\n'.join(rep),'\n'.join(gaprep),sep='\n')
    rgb=colors[g]
    lines=[]
    for nm,a,b,lab,sd,idn in bed:
        lines.append(f"{nm}\t{a}\t{b}\t{g} {lab} ({b-a} bp)\t0\t{sd}\t{a}\t{b}\t{rgb}")
    h1=seqs[0]; h2=seqs[1]
    j,a2,b2={'TMEM114':(15709,15708,79274),'SMCO1':(15904,15502,23296)}[g]; L=b2-a2
    lines.append(f"{names[0]}\t{max(j-1,0)}\t{j+1}\t{g} hap1 deletion junction (hap1 pos {j})\t0\t.\t{max(j-1,0)}\t{j+1}\t255,0,0")
    lines.append(f"{names[1]}\t{a2}\t{b2}\t{g} segment absent from hap1 ({b2-a2} bp)\t0\t.\t{a2}\t{b2}\t255,150,150")
    lines.append(f"{names[1]}\t{max(a2-1,0)}\t{a2+1}\t{g} breakpoint L\t0\t.\t{max(a2-1,0)}\t{a2+1}\t255,0,0")
    lines.append(f"{names[1]}\t{b2-1}\t{b2+1}\t{g} breakpoint R\t0\t.\t{b2-1}\t{b2+1}\t255,0,0")
    print(g,'DELETION: hap1 junction',j,'hap2 absent interval',a2,b2,b2-a2,'gapcols',L)
    open(D+g+'.bed','w').write('\n'.join(lines)+'\n')
    with open(D+g+'.aln.fa','w') as f:
        for nm,sq in zip(names,seqs):
            f.write('>'+nm+'\n'+'\n'.join(sq[i:i+80] for i in range(0,len(sq),80))+'\n')
    print('lens',[len(u) for u in ung])
