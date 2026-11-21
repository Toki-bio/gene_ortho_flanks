# Full locus: mouse 4600-6443 (left flank, the paper's 1228-bp segment 5001-6228, right flank);
# each dipodoid from its projection of mouse 4600 (as cut_leftend.py) to its right-end anchor + 300 (as cut_rightend.py).
# Writes full.fa and full.bed (per-row ungapped 0-based positions).
def rf(p): return "".join(l.strip().upper() for l in open(p,encoding="utf-8") if not l.startswith(">"))
left={"jaculus":(4833,7478),"zapus":(4817,7312),"bet_hap1":(4799,7280),"bet_hap2":(4799,7280),"tri_hap1":(4799,7280),"tri_hap2":(4799,7279)}
right={"jaculus":(6147,7710),"zapus":(6143,7656),"bet_hap1":(6145,7564),"bet_hap2":(6145,7565),"tri_hap1":(6145,7579),"tri_hap2":(6145,7569)}
anc={"jaculus":7702,"zapus":7656,"bet_hap1":7562,"bet_hap2":7563,"tri_hap1":7577,"tri_hap2":7567}
blk={"jaculus":308}  # block start within the right-end row (others 304), from rightend.bed
ren={"jaculus":"jerboa_Jjac","zapus":"Zapus_hud","bet_hap1":"Sbet_hap1","bet_hap2":"Sbet_hap2","tri_hap1":"Stri_hap1","tri_hap2":"Stri_hap2"}
m=rf("mouse_mm10_chr3_67261489-67272716.fa")
bed=[]
def B(n,s,e,lab,col): bed.append(f"{n}\t{s}\t{e}\t{lab}\t0\t+\t{s}\t{e}\t{col}")
with open("full.fa","w",newline="\n") as o:
    ms=m[4599:6443]; o.write(f">mouse_mm10\n{ms}\n")
    B("mouse_mm10",401,1629,"paper's 1228-bp segment (mm10 chr3:67266489-67267716)","230,160,0")
    B("mouse_mm10",401,1490,"deleted in all dipodoids (~1089 bp, mouse offsets 5001-6089)","200,0,0")
    B("mouse_mm10",1547,1570,"block kept in all (mouse 6147-6169)","0,120,110")
    assert ms[1547:1570]==m[6146:6169]
    for n in left:
        a=left[n][1]-(left[n][0]-4600); e=anc[n]+300
        j=right[n][1]-(right[n][0]-6090)          # 1-based window position of mouse 6090 = first base right of the junction
        s=rf(f"win/{n}.fa"); row=s[a-1:e]; o.write(f">{ren[n]}\n{row}\n")
        jb=j-a; bs=(anc[n]-301+blk.get(n,304))-(a-1)
        B(ren[n],max(jb-2,0),jb+2,"junction: deletion here (+-15 bp)","200,0,0")
        B(ren[n],bs,bs+23,"block kept in all (mouse 6147-6169)","0,120,110")
        print(n,len(row),"junction",jb,"block",row[bs:bs+23],"mouse",m[6146:6169])
open("full.bed","w",newline="\n").write("\n".join(bed)+"\n")
