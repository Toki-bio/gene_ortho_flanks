#!/bin/bash
# scfish: pair mode identical to the previous version; single-read mode + MIN_HITS on a known answer.
set -euo pipefail
D=/home/toki/gof/scfish
cd $D || exit 1
gcc -O3 -march=native -o scfish scfish.c
gcc -O3 -march=native -o scfish_old scfish_old.c
python3 - <<'PY'
import random
random.seed(7)
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))
anchors = [rnd(25) for _ in range(50)]
open("anchors.txt", "w").write("\n".join(anchors) + "\n")
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
with open("pairs.tsv", "w") as p, open("single.tsv", "w") as s, open("truth.txt", "w") as t:
    for i in range(2000):
        k = random.choice([0, 0, 0, 1, 2, 3])          # anchors planted in read i
        r = rnd(15000)
        for j in range(k):
            a = random.choice(anchors)
            a = rc(a) if random.random() < 0.5 else a
            pos = random.randint(0, 14900)
            r = r[:pos] + a + r[pos + 25:]
        s.write(f"@r{i}\t{r}\t+\t{'I'*len(r)}\n")
        # truth = anchor k-mer positions actually present (a later plant can overwrite an earlier one)
        aset = set(anchors) | {rc(a) for a in anchors}
        t.write(f"r{i}\t{sum(1 for x in range(len(r) - 24) if r[x:x + 25] in aset)}\n")
        r1, r2 = rnd(150), rnd(150)
        if i % 3 == 0:
            r2 = r2[:50] + random.choice(anchors) + r2[75:]
        p.write(f"@p{i}/1\t{r1}\t+\t{'I'*150}\t@p{i}/2\t{r2}\t+\t{'I'*150}\n")
PY
./scfish_old anchors.txt < pairs.tsv > old.out
./scfish anchors.txt < pairs.tsv > new.out
cmp old.out new.out && echo "pair mode identical ($(wc -l < new.out) pairs)"
for m in 1 2 3; do
  got=$(./scfish anchors.txt $m < single.tsv | cut -f1 | tr -d @ | sort)
  want=$(awk -v m=$m '$2>=m{print $1}' truth.txt | sort)
  [ "$got" = "$want" ] && echo "single MIN_HITS=$m: $(echo "$got" | wc -l) reads, matches truth" || echo "single MIN_HITS=$m: MISMATCH"
done
