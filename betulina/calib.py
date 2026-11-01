"""Shared helpers: read the gc_curve.py calibration, FASTA, GC."""
import statistics


def load_curve(path):
    """-> (expected(gc) -> relative two-copy depth, REF two-copy depth)"""
    ls = open(path).read().splitlines()
    exp = {int(g): float(v) for g, v in (x.split(":") for x in ls[0].split(":", 1)[1].split())}
    ref = float(ls[1].split("\t")[1])
    return (lambda g: exp[min(exp, key=lambda k: abs(k - g))]), ref


def read_fa(path):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def gc(x):
    n = sum(x.count(c) for c in "ACGT")
    return 100 * (x.count("G") + x.count("C")) / n if n else 45.0


def med(v):
    return statistics.median(v) if v else float("nan")
