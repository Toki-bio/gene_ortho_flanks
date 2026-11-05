# compare bash/awk rows with the Python controls.reads.tsv, read by read
# rows: id tag read call shared dt do cross cov x ind indbp big spans
# python: id read call shared dt do cross "own:cov=..,x=..,ind=..,indbp=..,big=..,spans=..;other:..."
FNR == NR {
	k = $1 SUBSEP $3
	call[k] = $4; sh[k] = $5; dt[k] = $6; dd[k] = $7; cr[k] = $8
	f = k SUBSEP $2; cov[f] = $9; x[f] = $10; ind[f] = $11; sp[f] = $14
	next
}
{
	k = $1 SUBSEP $2; n++
	if (!(k in call)) { miss++; next }
	if (call[k] != $3) { dc++; pair[$3 " -> " call[k]]++; if (shown++ < 12) print "call differs", $1, $2, "py", $3, $4, $5, $6, "awk", call[k], sh[k], dt[k], dd[k] }
	else sc++
	d = sh[k] - $4; if (d < 0) d = -d
	if (d > 100) dsh++
	if (dt[k] != $5 || dd[k] != $6) dev++
	if (cr[k] != $7) dcr++
	if ($7 == 1) {
		m = split($8, p, ";")
		for (i = 1; i <= m; i++) {
			split(p[i], q, ":"); tag = q[1]; split(q[2], v, ",")
			split(v[2], a, "="); split(v[6], b, "=")
			f = k SUBSEP tag; nf++
			if (x[f] != a[2]) { dx++; if (shown2++ < 8) print "fill x differs", $1, $2, tag, "py", a[2], "awk", x[f] }
			if (sp[f] != b[2]) dsp++
		}
	}
}
END {
	print "python reads", n, "missing in awk", miss + 0
	print "same call", sc + 0, "different call", dc + 0
	for (p2 in pair) print "  ", p2, pair[p2]
	print "shared differs >100 bp", dsh + 0, "| event counts differ", dev + 0, "| crossing differs", dcr + 0
	print "fill rows compared", nf + 0, "| substitutions differ", dx + 0, "| spans differs", dsp + 0
}
