import subprocess,glob
D='/home/toki/gof/views_het/'
MAFFT=glob.glob('/home/toki/miniforge3/envs/sinederella/bin/mafft')[0]
def readfa(p):
    n=[];s=[]
    for l in open(p):
        l=l.rstrip()
        if l.startswith('>'): n.append(l[1:]);s.append([])
        else: s[-1].append(l)
    return n,[''.join(x) for x in s]
def wr(p,n,s):
    with open(p,'w') as f:
        for a,b in zip(n,s): f.write('>'+a+'\n'+b+'\n')
n,s=readfa(D+'TMEM114.aln.raw.fa')
h1=s[0].replace('-','')
rows=s[1:]; names=n[1:]
# columns of hap2 breakpoints
def col_of(row,pos):  # column index (0-based) where ungapped count reaches pos
    c=0
    for i,ch in enumerate(row):
        if ch!='-':
            if c==pos: return i
            c+=1
    return len(row)
L,R=15708,79274
c1=col_of(rows[0],L); c2=col_of(rows[0],R)
print('cols',c1,c2)
h1L=h1[:15709]; h1R=h1[15709:]
res=[]
for tag,(a,b),piece in (('L',(0,c1),h1L),('M',(c1,c2),''),('R',(c2,len(rows[0])),h1R)):
    blk=[r[a:b] for r in rows]
    if not piece:
        res.append(['-'*(b-a)]+blk); continue
    # drop all-gap columns of the block for mafft input
    keep=[i for i in range(b-a) if any(r[i]!='-' for r in blk)]
    blk2=[''.join(r[i] for i in keep) for r in blk]
    wr(D+'blk.fa',names,blk2)
    open(D+'new.fa','w').write('>TMEM114_tri_hap1\n'+piece+'\n')
    out=subprocess.run([MAFFT,'--thread','8','--add',D+'new.fa','--quiet',D+'blk.fa'],capture_output=True,text=True,check=True).stdout
    on,os_=[],[]
    for l in out.split('\n'):
        if l.startswith('>'): on.append(l[1:]);os_.append([])
        elif l: os_[-1].append(l)
    d={a:''.join(b) for a,b in zip(on,os_)}
    res.append([d['TMEM114_tri_hap1']]+[d[x] for x in names])
final=[''.join(r[i] for r in res) for i in range(4)]
fn=['TMEM114_tri_hap1']+names
assert len(set(map(len,final)))==1
wr(D+'TMEM114.aln.raw.fa',fn,final)
print([len(x.replace('-','')) for x in final],len(final[0]))
