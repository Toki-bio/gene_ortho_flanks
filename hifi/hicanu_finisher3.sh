#!/bin/bash
# Finisher v3 (2026-11-04): tgStoreLoad -test showed tails of KILLED utgcns runs are
# corrupt (half-flushed last record). v2 bug: it appended the stalled attempt's tmp
# file to the part. v3 rules: (a) append attempt output ONLY on clean (rc=0, no stall)
# exit and move on to the next subrange; (b) pre-split ranges around the two known
# un-summable tigs 6925 / 5559 (both failed 8x rescue) so nothing re-stalls;
# (c) validate every part with tgStoreLoad -test before concat;
# (d) final cns from clean parts only: 0020 = a15+s+n, 0021 = a+s+n.
set -uo pipefail
cd /home/toki/gof/asm2/hicanu/unitigging/5-consensus || exit 1
export TMPDIR=/home/toki/tmp PATH=/home/toki/miniforge3/envs/canu/bin:$PATH
U=/home/toki/miniforge3/envs/canu/bin/utgcns
T=/home/toki/miniforge3/envs/canu/bin/tgStoreLoad
C="-T ../tri.ctgStore 1 -e 0.05 -em 0.05 -pbdagcon -edlib -threads 8 -maxcoverage 15"
echo "#### finisher3 start $(date '+%F %T')"
rm -f ctgcns/0020n.cns ctgcns/0021n.cns ctgcns/0020.cns ctgcns/0021.cns

watch_ranges() {  # part outbase range...
  PART=$1; OB=$2; shift 2
  for R in "$@"; do
    att=0
    while [ $att -lt 6 ]; do
      att=$((att+1))
      L=$OB.a$att.log; MT=$OB.a$att.tmp.cns
      rm -f ctgcns/$MT
      timeout 7200 stdbuf -oL $U -R ../tri.ctgStore/partition.$PART -P $PART $C \
        -tig "$R" -O ctgcns/$MT > $L 2>&1 &
      PID=$!; s1=0; s2=0; s3=0; s4=0; stalled=0
      while kill -0 $PID 2>/dev/null; do
        sleep 300
        s4=$s3; s3=$s2; s2=$s1; s1=$(stat -c%s ctgcns/$MT 2>/dev/null || echo 0)
        if [ "$s1" = "$s2" ] && [ "$s2" = "$s3" ] && [ "$s3" = "$s4" ]; then
          stalled=1; echo "  $OB range $R attempt $att STALL (size $s1)"; kill $PID; break
        fi
      done
      wait $PID; rc=$?
      if [ $stalled -eq 0 ] && [ $rc -eq 0 ]; then
        cat ctgcns/$MT >> ctgcns/$OB.cns; rm -f ctgcns/$MT
        echo "  $OB range $R OK $(date '+%T')"; break
      fi
      rm -f ctgcns/$MT     # stalled attempt output is poison: never append it
      STIG=$(tail -1 $L | sed 's/^ *//' | cut -d' ' -f1)
      case $STIG in (*[!0-9]*|"") echo "  $OB range $R: cannot parse stalled tig, ABORT"; return 1;; esac
      echo "  $OB range $R stalled on tig $STIG (dropped), retrying from $((STIG+1))"
      R="$((STIG+1))-99999"
    done
  done
  return 0
}

( watch_ranges 0020 0020n 6470-6924 6926-99999; echo $? > 0020n.rc ) &
W20=$!
( watch_ranges 0021 0021n 5110-5558 5560-99999; echo $? > 0021n.rc ) &
W21=$!
wait $W20; echo "p20 ranges rc=$(cat 0020n.rc) $(date '+%T')"
wait $W21; echo "p21 ranges rc=$(cat 0021n.rc) $(date '+%T')"

echo "#### validating all parts with tgStoreLoad -test"
ok=1
for f in ctgcns/0020a15.cns ctgcns/0020s.cns ctgcns/0020n.cns ctgcns/0021a.cns ctgcns/0021s.cns ctgcns/0021n.cns; do
  $T -S ../tri.seqStore -T ../tri.ctgStore 2 -test "$f" >/dev/null 2>$f.test.err
  if [ $? -eq 0 ]; then echo "  $f CLEAN"; else echo "  $f CORRUPT"; ok=0; fi
done
[ $ok -eq 1 ] || { echo "ABORT: corrupt part present, not building cns"; exit 1; }

cat ctgcns/0020a15.cns ctgcns/0020s.cns ctgcns/0020n.cns > ctgcns/0020.cns
cat ctgcns/0021a.cns   ctgcns/0021s.cns ctgcns/0021n.cns > ctgcns/0021.cns
echo "#### built 0020.cns $(stat -c%s ctgcns/0020.cns) B; 0021.cns $(stat -c%s ctgcns/0021.cns) B"
# sanity: also test the two assembled partition files themselves
for p in 0020 0021; do
  $T -S ../tri.seqStore -T ../tri.ctgStore 2 -test ctgcns/$p.cns >/dev/null 2>&1 && echo "  $p.cns CLEAN" || { echo "  $p.cns CORRUPT"; exit 1; }
done

echo "#### resuming canu $(date '+%F %T')"
cd /home/toki/gof/asm2 || exit 1
canu -p tri -d hicanu genomeSize=2.8g useGrid=false maxThreads=56 maxMemory=240g \
  -pacbio-hifi /home/toki/gof/reads/ERR16751478.fastq.gz > hicanu_resume.out 2> hicanu_resume.err
echo "#### canu resume exit $? $(date '+%F %T')"
if [ -s hicanu/tri.contigs.fasta ]; then
  awk '/^>/{n++} !/^>/{s+=length($0)} END{print n" contigs "s" bp"}' hicanu/tri.contigs.fasta
  /home/toki/miniforge3/envs/asm25/bin/gfastats hicanu/tri.contigs.fasta 2>/dev/null | head -n 20
else
  echo "#### ERROR: no tri.contigs.fasta"
  tail -5 hicanu_resume.err
fi
df -h /home | tail -1
echo "#### finisher3 done $(date '+%F %T')"
