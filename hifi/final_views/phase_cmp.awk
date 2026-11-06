# Phase test of one fill. Input: three PAF files with cs tags, reads mapped alone to T (target window, N kept), O (other haplotype), F (fixed window).
# usage: awk -F'\t' -v g0= -v g1= -v f0= -v f1= -v id= -v tag= -v Tf=T.paf -v Of=O.paf -v Ff=F.paf -v rows=rows.tsv -f phase_cmp.awk T.paf O.paf F.paf
# Read assigned by difference events over the same read bases in the flanks, then substitutions inside the fill.
# prints: reads T_reads O_reads tie foreign short one_side T_cross O_cross tie_cross T_span T_x T_x_median T_indel_median T_big
#         O_x_median tie_x_median T_enter T_enter_x_per_kb O_enter_x_per_kb verdict
BEGIN { EDGE = 500; EMIN = 400; BIG = 50; SHARED = 2000; MARGIN = 3; OWN = 0.005; FAR = 1000000000000 }

function min(a, b) { return a < b ? a : b }
function max(a, b) { return a > b ? a : b }
function getcs(   i) { for (i = 13; i <= NF; i++) if (substr($i, 1, 5) == "cs:Z:") return substr($i, 6); return "" }

# matched bases r..r+n-1 at low read coordinate q: update read span left and right of the hole, and edge counts
function seg(r, n, q, h0, h1, fw, st, rd,    lf, s, off, a, b) {
	lf = min(r + n, h0) - r
	if (lf > 0) {
		if (fw) { a = q; b = q + lf } else { a = q + n - lf; b = q + n }
		if (llo < 0 || a < llo) llo = a
		if (b > lhi) lhi = b
		LE[rd, st] += max(0, min(r + n, h0) - max(r, h0 - EDGE))
	}
	s = max(r, h1)
	if (r + n > s) {
		off = s - r
		if (fw) { a = q + off; b = q + n } else { a = q; b = q + n - off }
		if (rlo < 0 || a < rlo) rlo = a
		if (b > rhi) rhi = b
		RE[rd, st] += max(0, min(r + n, h1 + EDGE) - max(r, h1))
	}
}

# close the read stretch being built: inserted bases and N positions are not aligned bases
function flush(INT, rd) {
	if (lhi > llo) INT[rd] = INT[rd] " " llo ":" lhi
	if (rhi > rlo) INT[rd] = INT[rd] " " rlo ":" rhi
	llo = -1; lhi = -1; rlo = -1; rhi = -1
}

# walk one alignment on T or O: read spans into INT[rd], events into EV[rd]
function flank(h0, h1, INT, EV,    rd, st, qs, qe, r, o, fw, cs, c, n, q, rb, qb) {
	rd = $1; st = $5; qs = $3 + 0; qe = $4 + 0; r = $8 + 0; o = 0; fw = (st == "+"); cs = getcs()
	llo = -1; lhi = -1; rlo = -1; rhi = -1
	while (cs != "") {
		c = substr(cs, 1, 1)
		if (c == ":") {
			match(cs, /^:[0-9]+/); n = substr(cs, 2, RLENGTH - 1) + 0; cs = substr(cs, RLENGTH + 1)
			seg(r, n, fw ? qs + o : qe - o - n, h0, h1, fw, st, rd)
			r += n; o += n
		} else if (c == "*") {
			rb = substr(cs, 2, 1); qb = substr(cs, 3, 1); cs = substr(cs, 4)
			q = fw ? qs + o : qe - o - 1
			if (rb != "n" && qb != "n") {
				seg(r, 1, q, h0, h1, fw, st, rd)
				if (!(r >= h0 && r < h1)) EV[rd] = EV[rd] " " q
			} else flush(INT, rd)
			r++; o++
		} else {
			match(cs, /^[-+][a-z]+/); n = RLENGTH - 1; cs = substr(cs, RLENGTH + 1)
			if (c == "+") q = fw ? qs + o : qe - o - n; else q = fw ? qs + o : qe - o
			if (n < BIG && !(r >= h0 - 1 && r <= h1)) EV[rd] = EV[rd] " " max(q - 1, 0)
			if (c == "+") { flush(INT, rd); o += n } else r += n
		}
	}
	flush(INT, rd)
}

function fillscan(    rd, qs, r, o, cs, c, n, rb, qb, re_) {
	rd = $1; r = $8 + 0; re_ = $9 + 0; cs = getcs()
	if (r <= f0 - EMIN && re_ >= f1 + EMIN) SP[rd] = 1
	while (cs != "") {
		c = substr(cs, 1, 1)
		if (c == ":") {
			match(cs, /^:[0-9]+/); n = substr(cs, 2, RLENGTH - 1) + 0; cs = substr(cs, RLENGTH + 1)
			FC[rd] += max(0, min(r + n, f1) - max(r, f0)); r += n
		} else if (c == "*") {
			rb = substr(cs, 2, 1); qb = substr(cs, 3, 1); cs = substr(cs, 4)
			if (rb != "n" && qb != "n" && r >= f0 && r < f1) { FX[rd]++; FC[rd]++ }
			r++
		} else {
			match(cs, /^[-+][a-z]+/); n = RLENGTH - 1; cs = substr(cs, RLENGTH + 1)
			if (r >= f0 && r < f1) {
				FI[rd]++; FB[rd] += n
				if (n >= BIG) FG[rd] = FG[rd] (FG[rd] == "" ? "" : "/") c n
			}
			if (c == "-") r += n
		}
	}
}

# sorted, merged spans from "a:b a:b" into A[1..],B[1..]; returns count
function spans(s, A, B,    k, n, i, j, p, ta, tb, m) {
	n = split(s, p, " ")
	for (i = 1; i <= n; i++) { split(p[i], m, ":"); ta[i] = m[1] + 0; tb[i] = m[2] + 0 }
	for (i = 2; i <= n; i++) {
		a = ta[i]; b = tb[i]
		for (j = i - 1; j >= 1 && ta[j] > a; j--) { ta[j + 1] = ta[j]; tb[j + 1] = tb[j] }
		ta[j + 1] = a; tb[j + 1] = b
	}
	k = 0
	for (i = 1; i <= n; i++) {
		if (k && ta[i] <= B[k]) B[k] = max(B[k], tb[i])
		else { k++; A[k] = ta[i]; B[k] = tb[i] }
	}
	return k
}

function inside(p, A, B, n,    i) { for (i = 1; i <= n; i++) if (p >= A[i] && p < B[i]) return 1; return 0 }

function count(s, A, B, n, C, D, m,    k, i, e, c) {
	k = split(s, e, " "); c = 0
	for (i = 1; i <= k; i++) if (inside(e[i] + 0, A, B, n) && inside(e[i] + 0, C, D, m)) c++
	return c
}

function median(s,    n, v, i, j, a) {
	n = split(s, v, " ")
	if (!n) return "NA"
	for (i = 1; i <= n; i++) v[i] += 0
	for (i = 2; i <= n; i++) {
		a = v[i]
		for (j = i - 1; j >= 1 && v[j] > a; j--) v[j + 1] = v[j]
		v[j + 1] = a
	}
	return n % 2 ? v[(n + 1) / 2] : (v[n / 2] + v[n / 2 + 1]) / 2
}

FILENAME == Tf { flank(g0, g1, TI, TE); seen[$1] = 1; next }
FILENAME == Of { flank(FAR, FAR, OI, OE); next }
FILENAME == Ff { fillscan(); next }

END {
	for (rd in seen) {
		delete A; delete B; delete C; delete D
		n = spans(TI[rd], A, B); shared = 0; dt = 0; dd = 0
		if (!(rd in OI)) call = "one_side"
		else {
			m = spans(OI[rd], C, D)
			for (i = 1; i <= n; i++) for (j = 1; j <= m; j++) shared += max(0, min(B[i], D[j]) - max(A[i], C[j]))
			dt = count(TE[rd], A, B, n, C, D, m); dd = count(OE[rd], A, B, n, C, D, m)
			if (shared < SHARED) call = "short"
			else if (min(dt, dd) / shared > OWN) call = "foreign"
			else if (dd - dt >= MARGIN && dt * 2 <= dd) call = "T"
			else if (dt - dd >= MARGIN && dd * 2 <= dt) call = "O"
			else call = "tie"
		}
		calls[call]++; tot++
		st = (LE[rd, "+"] + RE[rd, "+"] >= LE[rd, "-"] + RE[rd, "-"]) ? "+" : "-"
		x = (LE[rd, st] >= EMIN && RE[rd, st] >= EMIN) ? 1 : 0
		cov = FC[rd] + 0; xs = FX[rd] + 0; ind = FI[rd] + 0; bp = FB[rd] + 0; sp = SP[rd] + 0
		if (x) {
			cross[call]++
			if (sp) {
				XS[call] = XS[call] " " xs; IS[call] = IS[call] " " ind; NS[call]++
				if (FG[rd] != "") BG[call] = BG[call] (BG[call] == "" ? "" : "/") FG[rd]
			}
		} else if (cov >= EDGE) { en[call]++; ec[call] += cov; ex[call] += xs }
		print id, tag, rd, call, shared, dt, dd, x, cov, xs, ind, bp, FG[rd], sp > rows
	}
	if (cross["T"] + 0 == 0) {
		if (cross["O"] + 0 == 0 && cross["tie"] + 0 > 0) v = "flanks_tie"
		else if (calls["T"] + calls["O"] == 0) v = calls["tie"] + 0 ? "flanks_tie" : "no_phasable_reads"
		else v = "no_target_cross"
	} else if (NS["T"] + 0 == 0) v = "wrong_allele_structure"
	else { md = median(XS["T"]); v = md == 0 ? "own_or_identical" : (md >= 2 ? "wrong_allele" : "undecided") }
	tx = XS["T"]; gsub(/^ /, "", tx); gsub(/ /, ",", tx)
	te = ec["T"] ? sprintf("%.2f", 1000 * ex["T"] / ec["T"]) : "NA"
	oe = ec["O"] ? sprintf("%.2f", 1000 * ex["O"] / ec["O"]) : "NA"
	print tot + 0, calls["T"] + 0, calls["O"] + 0, calls["tie"] + 0, calls["foreign"] + 0, calls["short"] + 0, calls["one_side"] + 0, \
		cross["T"] + 0, cross["O"] + 0, cross["tie"] + 0, NS["T"] + 0, tx, median(XS["T"]), median(IS["T"]), BG["T"], \
		median(XS["O"]), median(XS["tie"]), en["T"] + 0, te, oe, v
}
