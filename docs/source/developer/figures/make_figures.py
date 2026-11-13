#!/usr/bin/env python3
"""Draw the class scheme of the Developer Guide from its two tables.

    python3 docs/source/developer/figures/make_figures.py

reads classes.tsv (the boxes) and arrows.tsv (the arrows) beside this file and writes, for every figure named in
them, <figure>.dot (the Graphviz source) and <figure>.svg (drawn by Graphviz's dot: the program named by the
environment variable DOT, else the first dot on PATH).  The site's build uses the committed SVGs; it needs no Graphviz.

A change of structure is an edit of the two tables and a run of this script; look at the SVGs afterwards.  The
consistency test (tests/consistency/repo_consistency_test.py, check 44) checks that every type of classes.tsv is
declared in the file the table names.  Standard library only.
"""
import csv
import html
import os
import shutil
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))

TITLES = {
    "libraries": "Libraries and their link lines",
    "setup": "Project and set-up",
    "evaluation": "Evaluation",
    "inference": "Inference: objective, parameters, fitters, sampler",
    "frontends": "Front ends",
}

# Box styles by state: immutable after the set-up (double border), written by every evaluation (rounded, grey),
# neither (plain).  Legible in black and white.
NODE_STYLE = {
    "immutable": 'shape=box, peripheries=2, style="filled", fillcolor="white"',
    "evaluation": 'shape=box, style="rounded,filled", fillcolor="#d9d9d9"',
    "run": 'shape=box, style="filled", fillcolor="white"',
}
# Arrow styles by kind.
EDGE_STYLE = {
    "owns": 'style=solid, dir=both, arrowtail=diamond, arrowhead=normal',
    "reads": 'style=dashed, arrowhead=vee',
    "makes": 'style=dotted, penwidth=1.6, arrowhead=empty',
    "links": 'style=solid, arrowhead=normal',
}
FONT = 'fontname="Helvetica"'
JOB_WIDTH = 22  # characters per line of a box's job
# The libraries flow top to bottom (a library above the ones it links); a layer left to right (owner left of owned),
# which keeps the figures narrow enough for the page.
RANKDIR = {"libraries": "TB"}


def read_table(name):
    with open(os.path.join(HERE, name), newline="") as f:
        rows = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def node_id(name):
    return "n_" + "".join(c if c.isalnum() else "_" for c in name)


def node_label(row):
    name = html.escape(row["type"])
    job = "<br/>".join(html.escape(line) for line in textwrap.wrap(row["job"], JOB_WIDTH))
    return f'<<b>{name}</b><br/><font point-size="12">{job}</font>>'


def legend(figure):
    """The legend: one node, a table of the box and arrow kinds (a node, so that it takes no width of the layout)."""
    if figure == "libraries":
        rows = '<TR><TD ALIGN="LEFT">an arrow from A to B: A links B (CMake) or imports it (Python)</TD></TR>'
        return [f'  legend [shape=plaintext, fontsize=11, label=<<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="2" '
                f'COLOR="#888888"><TR><TD><B>Legend</B></TD></TR>{rows}</TABLE>>];']
    box = '<TABLE BORDER="1" CELLBORDER="0" CELLPADDING="3" {extra}><TR><TD>{text}</TD></TR></TABLE>'
    immutable = ('<TABLE BORDER="1" CELLBORDER="1" CELLSPACING="2" CELLPADDING="2"><TR><TD>immutable after the set-up'
                 '</TD></TR></TABLE>')
    evaluation = box.format(extra='STYLE="rounded" BGCOLOR="#d9d9d9"', text="written by every evaluation")
    run = box.format(extra="", text="a run or session object")
    arrows = [("owns", "solid line, diamond at the owner"),
              ("reads", "dashed line; label nc: not const"),
              ("makes", "dotted line, open head")]
    rows = "".join(f'<TR><TD ALIGN="LEFT"><B>{k}</B></TD><TD ALIGN="LEFT">{t}</TD></TR>' for k, t in arrows)
    table = ('<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="3" COLOR="#888888">'
             '<TR><TD COLSPAN="2"><B>Legend</B></TD></TR>'
             f'<TR><TD COLSPAN="2">{immutable}</TD></TR><TR><TD COLSPAN="2">{evaluation}</TD></TR>'
             f'<TR><TD COLSPAN="2">{run}</TD></TR>{rows}</TABLE>')
    return [f"  legend [shape=plaintext, fontsize=11, label=<{table}>];"]


def dot_source(figure, boxes, arrows):
    rankdir = RANKDIR.get(figure, "LR")
    lines = [f"digraph {figure} {{",
             f'  graph [rankdir={rankdir}, {FONT}, fontsize=16, labelloc=t, label="{TITLES.get(figure, figure)}", '
             'nodesep=0.18, ranksep=0.3, newrank=true];',
             f"  node [{FONT}, fontsize=14, margin=\"0.08,0.04\"];",
             f"  edge [{FONT}, fontsize=12];"]
    for row in boxes:
        lines.append(f'  {node_id(row["type"])} [label={node_label(row)}, {NODE_STYLE[row["state"]]}];')
    for a in arrows:
        style = EDGE_STYLE[a["kind"]]
        extra = ', label=" nc"' if a["kind"] == "reads" and a["access"] == "nonconst" else ""
        lines.append(f'  {node_id(a["from"])} -> {node_id(a["to"])} [{style}{extra}];')
    lines += legend(figure)
    lines.append("}")
    return "\n".join(lines) + "\n"


def main():
    classes = read_table("classes.tsv")
    arrows = read_table("arrows.tsv")
    problems = []
    for r in classes:
        if r["state"] not in NODE_STYLE:
            problems.append(f'classes.tsv: {r["type"]} ({r["figure"]}): unknown state {r["state"]!r}')
    figures = []
    for r in classes:
        if r["figure"] not in figures:
            figures.append(r["figure"])
    for a in arrows:
        if a["kind"] not in EDGE_STYLE:
            problems.append(f'arrows.tsv: {a["from"]} -> {a["to"]}: unknown kind {a["kind"]!r}')
        if (a["kind"] == "links") != (a["figure"] == "libraries"):
            problems.append(f'arrows.tsv: {a["from"]} -> {a["to"]}: kind links belongs to the libraries figure only')
        if a["kind"] == "reads" and a["access"] not in ("const", "nonconst"):
            problems.append(f'arrows.tsv: {a["from"]} -> {a["to"]}: a reads arrow needs access const or nonconst')
        if not a["where"].strip():
            problems.append(f'arrows.tsv: {a["from"]} -> {a["to"]}: no "where"')
        names = {r["type"] for r in classes if r["figure"] == a["figure"]}
        for end in (a["from"], a["to"]):
            if end not in names:
                problems.append(f'arrows.tsv: {end} is no box of the figure {a["figure"]} in classes.tsv')
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    dot = os.environ.get("DOT") or shutil.which("dot")
    if not dot:
        print("make_figures.py: Graphviz's dot not found (set DOT)", file=sys.stderr)
        return 1
    for figure in figures:
        boxes = [r for r in classes if r["figure"] == figure]
        edges = [a for a in arrows if a["figure"] == figure]
        source = dot_source(figure, boxes, edges)
        base = os.path.join(HERE, figure)
        with open(base + ".dot", "w") as f:
            f.write(source)
        svg = subprocess.run([dot, "-Tsvg", base + ".dot"], check=True, capture_output=True, text=True).stdout
        with open(base + ".svg", "w") as f:
            f.write(svg)
        print(f"{figure}: {len(boxes)} boxes, {len(edges)} arrows -> {os.path.relpath(base, HERE)}.svg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
