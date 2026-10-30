#!/usr/bin/env python3
"""Positions of N runs in a query interval and distance to a site. Usage: nrun_pos.py Q2BIT CHR:S-E SITE_POS"""
import re, subprocess, sys
q2bit, loc, site = sys.argv[1], sys.argv[2], int(sys.argv[3])
c, se = loc.split(":"); s, e = map(int, se.split("-"))
seq = "".join(subprocess.run(["twoBitToFa", f"{q2bit}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True).stdout.splitlines()[1:]).upper()
for m in re.finditer("N+", seq):
    a = s + m.start()
    print(f"N run {c}:{a}-{a + len(m.group())} len {len(m.group())}, distance to site {abs(a - site)} bp")
