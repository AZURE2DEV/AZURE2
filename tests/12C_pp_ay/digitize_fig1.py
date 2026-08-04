#!/usr/bin/env python3
"""Regenerate data/baumann_ay.dat from fig. 1 of Baumann et al. (1992).

    python3 digitize_fig1.py "/path/to/Baumann et al. (1992).pdf"

Renders page 3 at 600 dpi with pdftoppm, calibrates the axes on the printed
tick labels, and traces the drawn curve in each of the six panels. See
README.md for what the resulting numbers are and are not.
"""
import json, subprocess, sys, tempfile, os
import numpy as np
from PIL import Image
from scipy import ndimage

# Panel frames, located by looking for long horizontal and vertical dark runs
# in the 600 dpi render: (x0, x1, y0, y1, E_lab / MeV).
PANELS = [(1247, 2464,  679, 1586, 1.618), (2543, 3753,  679, 1586, 1.658),
          (1247, 2464, 1686, 2593, 1.708), (2543, 3753, 1686, 2593, 1.738),
          (1247, 2464, 2691, 3595, 1.758), (2543, 3753, 2691, 3595, 1.779)]
M = 26                                # margin clearing the interior tick marks
DEG_PER_PX = 150.0 / (1031 - 7)       # x labels 0 and 150 degrees
X_ZERO = 7.0
Y_ZERO = 449.0                        # the dashed A_y = 0 rule
PX_PER_UNIT = 93.125 / 0.2            # y ticks are 0.2 apart
ANGLES = list(range(40, 161, 10))
def runs_of(col):
    idx=np.where(col)[0]
    if idx.size==0: return []
    out=[]; s=idx[0]; p=idx[0]
    for v in idx[1:]:
        if v-p>2: out.append((s,p)); s=v
        p=v
    out.append((s,p)); return out

def trace(x0,x1,y0,y1):
    P=dark[y0+M:y1-M, x0+M:x1-M]
    # kill the thin dashed zero line (and any thin horizontal rule); the curve
    # is ~12 px thick so it survives a 7 px vertical opening
    P=ndimage.binary_opening(P, structure=np.ones((7,1),bool))
    lab,n=ndimage.label(P,structure=np.ones((3,3)))
    # the curve can fragment; keep every elongated piece and drop the compact
    # blobs, which are the "1.618 MeV" label characters
    C=np.zeros_like(P)
    for sl,i in zip(ndimage.find_objects(lab), range(1,n+1)):
        h=sl[0].stop-sl[0].start; w=sl[1].stop-sl[1].start
        if w>150 or h>150: C |= (lab==i)
    W=C.shape[1]; ymid=np.full(W,np.nan)
    # seed at the column with a single unambiguous run, then follow continuity
    seed=None
    for j in range(W):
        r=runs_of(C[:,j])
        if len(r)==1: seed=j; break
    if seed is None: return ymid
    def walk(rng, prev):
        for j in rng:
            r=runs_of(C[:,j])
            if not r: continue
            cs=[0.5*(s+e) for s,e in r]
            k=int(np.argmin([abs(c-prev) for c in cs]))
            if abs(cs[k]-prev) > 170:         # implausible jump: keep previous
                continue
            ymid[j]=cs[k]; prev=cs[k]
        return prev
    r0=runs_of(C[:,seed]); ymid[seed]=0.5*(r0[0][0]+r0[0][1])
    walk(range(seed+1,W), ymid[seed])
    walk(range(seed-1,-1,-1), ymid[seed])
    return ymid


def main(pdf):
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdftoppm", "-r", "600", "-f", "3", "-l", "3", "-png",
                        pdf, os.path.join(td, "p")], check=True)
        png = [f for f in os.listdir(td) if f.endswith(".png")][0]
        global dark
        dark = np.asarray(Image.open(os.path.join(td, png)).convert("L")) < 128

    rows = []
    for (x0, x1, y0, y1, E) in PANELS:
        ymid = trace(x0, x1, y0, y1)
        for th in ANGLES:
            j = int(round(X_ZERO + th / DEG_PER_PX - M))
            w = [ymid[k] for k in range(max(0, j - 3), min(len(ymid), j + 4))
                 if not np.isnan(ymid[k])]
            if not w:
                continue                      # steep section the trace lost
            ay = -((np.mean(w) + M) - Y_ZERO) / PX_PER_UNIT
            rows.append((E, float(th), round(float(ay), 3), 0.04))

    with open("data/baumann_ay.dat", "w") as f:
        for E, th, v, e in rows:
            f.write("%.4fe+00 %5.1f % .4f %.4f\n" % (E, th, v, e))
    print("wrote %d points" % len(rows))

if __name__ == "__main__":
    main(sys.argv[1])
