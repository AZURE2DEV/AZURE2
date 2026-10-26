#!/usr/bin/env python3
"""Write 13N_fast.azr: tests/13N/13N.azr without its capture segments (exit pair 2: the external
capture integrals are what make a 13N evaluation cost 0.7 s), every <targetInt> line inactive,
and the normalization of the meyer_84 elastic segment free with a 5 % systematic error
(segment-line fields 10 and 11): 314 elastic points, an evaluation in milliseconds, one
normalization for the chain to sample.  Regenerate with  python3 make_model.py  from here."""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
text = open(os.path.join(HERE, "..", "13N", "13N.azr")).read()

def edit_block(text, tag, fn):
    a = text.index("<%s>" % tag) + len("<%s>" % tag)
    b = text.index("</%s>" % tag)
    body = "\n".join(l for l in (fn(line) for line in text[a:b].split("\n")) if l is not None)
    return text[:a] + body + text[b:]

def seg_line(line):
    t = line.split()
    if len(t) >= 15 and t[2] == "2": return None      # a capture segment: dropped
    if len(t) >= 15 and t[14].endswith("meyer_84.dat"):
        t[9] = "1"; t[10] = "5"   # vary norm, 5 % systematic error
        return " ".join(t)
    return line

def ti_line(line):
    t = line.split()
    if len(t) > 2 and t[0] == "1": t[0] = "0"
    return " ".join(t)

text = edit_block(text, "segmentsData", seg_line)
text = edit_block(text, "targetInt", ti_line)
out = os.path.join(HERE, "13N_fast.azr")
open(out, "w").write(text)
print("wrote", out)
