#!/usr/bin/env python3
"""Independent chi2 of a covariance run of tests/data_covariance.

    python3 verify.py RUNDIR

RUNDIR holds t.azr, the covariance file, data/ and the output/ of a "calculate with data" run.
This script parses the covariance file itself, picks the fitted rows from the segment lines
(energy and angle windows), takes the residuals r_i = (fit - n d)/(n sigma) of those points
from output/AZUREOut_aa=1_R=1.out, and prints

    chi2_python  chi2_engine(covariance.out)  max|segment partition - chiSquared.out| (relative)

chi2_python = v^T C^-1 v with v_i = sigma_i r_i in the data's units (sigma_i = sqrt(C_ii) with the
data in place of the model in relvector terms, as the engine sets it) and, for relvector terms,
the model frozen in output/covarianceModel.out.
"""
import os
import sys

import numpy as np


def records(path):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t or t[0].startswith("#"):
                continue
            rows.append([float(x) for x in t[:4]])
    return rows


def read_cov(path):
    toks = []
    with open(path) as f:
        for line in f:
            toks += line.split("#", 1)[0].split()
    files, k = [], 0
    while toks[k] == "file":
        files.append((toks[k + 1], int(toks[k + 2])))
        k += 3
    N = sum(c for _, c in files)
    C = np.zeros((N, N))
    rel = []
    while k < len(toks):
        kind = toks[k]
        k += 1
        if kind == "block":
            n = int(toks[k]); k += 1
            idx = np.array([int(x) - 1 for x in toks[k:k + n]]); k += n
            B = np.array([float(x) for x in toks[k:k + n * n]]).reshape(n, n); k += n * n
            C[np.ix_(idx, idx)] += 0.5 * (B + B.T)
        elif kind == "diagonal":
            C[np.diag_indices(N)] += np.array([float(x) for x in toks[k:k + N]]); k += N
        elif kind == "vector":
            v = np.array([float(x) for x in toks[k:k + N]]); k += N
            C += np.outer(v, v)
        elif kind == "relvector":
            rel.append(np.array([float(x) for x in toks[k:k + N]])); k += N
        elif kind == "matrix":
            B = np.array([float(x) for x in toks[k:k + N * N]]).reshape(N, N); k += N * N
            C += 0.5 * (B + B.T)
        else:
            raise SystemExit("unknown term " + kind)
    return files, C, rel


def output_blocks(path):
    blocks, cur = [], []
    with open(path) as f:
        for line in f:
            t = line.split()
            if t:
                cur.append([float(x) for x in t])
            elif cur:
                blocks.append(cur)
                cur = []
    if cur:
        blocks.append(cur)
    return blocks


def main(run):
    lines = []
    text = open(os.path.join(run, "t.azr")).read()
    seg = text[text.index("<segmentsData>") + 15:text.index("</segmentsData>")]
    for line in seg.splitlines():
        t = line.split()
        if t and t[0] == "1":
            lines.append(t)
    covname = next(t[t.index("covariance") + 1] for t in lines if "covariance" in t)
    files, C, rel = read_cov(os.path.join(run, covname))
    offset, o = {}, 0
    for name, c in files:
        offset[name] = o
        o += c
    blocks = output_blocks(os.path.join(run, "output", "AZUREOut_aa=1_R=1.out"))
    rows, r, segof = [], [], []
    b = 0
    for s, t in enumerate(lines, start=1):
        if (t[1], t[2]) != ("1", "1"):
            continue
        block = blocks[b]
        b += 1
        if "covariance" not in t:
            continue
        e0, e1, a0, a1 = map(float, t[3:7])
        path = t[14]
        recs = records(os.path.join(run, path))
        sel = [k for k, x in enumerate(recs) if e0 <= x[0] <= e1 and a0 <= x[1] <= a1]
        if len(sel) != len(block):
            raise SystemExit("segment %d: %d records selected, %d in the output" % (s, len(sel), len(block)))
        for k, out in zip(sel, block):
            rows.append(offset[path] + k)
            r.append((out[3] - out[5]) / out[6])  # fit, data * n, error * n
            segof.append(s)
    rows = np.array(rows)
    r = np.array(r)
    data = np.zeros(C.shape[0])
    for name, c in files:
        data[offset[name]:offset[name] + c] = [x[2] for x in records(os.path.join(run, name))]
    sigma2 = np.diag(C).copy()
    for f in rel:
        sigma2 += (f * data) ** 2
    sigma = np.sqrt(sigma2[rows])
    Cs = C[np.ix_(rows, rows)].copy()
    if rel:
        m = np.full(C.shape[0], np.nan)
        with open(os.path.join(run, "output", "covarianceModel.out")) as f:
            for line in f:
                t = line.split()
                if t and not t[0].startswith("#"):
                    m[int(t[0]) - 1] = float(t[5])
        ms = m[rows]
        if np.isnan(ms).any():
            raise SystemExit("covarianceModel.out lacks fitted rows")
        for fr in rel:
            u = fr[rows] * ms
            Cs += np.outer(u, u)
    v = sigma * r
    z = np.linalg.solve(Cs, v)
    chi2 = float(v @ z)
    engine = None
    with open(os.path.join(run, "output", "covariance.out")) as f:
        for line in f:
            if not line.startswith("#"):
                engine = float(line.split(",")[4])
    # per-segment partition sum r_i (M^-1 r)_i = v_i z_i against chiSquared.out (6 digits)
    part = {}
    for s, x in zip(segof, v * z):
        part[s] = part.get(s, 0.0) + x
    worst = 0.0
    with open(os.path.join(run, "output", "chiSquared.out")) as f:
        for line in f:
            t = line.strip().split(",")
            if len(t) >= 2 and t[0].isdigit() and int(t[0]) in part:
                p = part[int(t[0])]
                worst = max(worst, abs(float(t[1]) - p) / max(abs(p), 1.0))
    print("%.12e %.12e %.3e" % (chi2, engine, worst))


if __name__ == "__main__":
    main(sys.argv[1])
