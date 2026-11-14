#!/usr/bin/env python3
"""A stand-in for the ENDF checking codes in the tests of pyazr.endf (not a checker: it checks nothing).

Linked under a code's name (checkr, fizcon, psyche, stanef, fixup, dictin, endf2c), it reads that code's deck as the
real program does -- on standard input the input file, the message (or output) file, the options, then the next
input file or DONE (checkr.f lines 255-270, 806-861; stanef.f lines 145-180); PREPRO's codes ENDFB.IN, ENDFB.OUT and
<CODE>.INP in the working directory (dictin.f lines 205-256) -- and writes what the real program writes, in its
layout, for the role named by $ENDF_STANDIN_ROLE:

- clean:    a report with no finding (CHECKR's "No problems to report", as checkr.out-0 of its sources);
- findings: a report with two findings and the closing count, lines as the built CHECKR 8.24 wrote them for the
            ENDF/B-VIII.1 3He + 4He file (2026-11-14), and PSYCHE's "FILE / SECTION" form;
- crash:    half a report, then a Fortran runtime error on standard error and exit status 2;
- lossy:    (rewriting codes) the file with every number of MF=3 rounded to 7 significant digits in the classic
            form, as STANEF's standard form does, and the sequence numbers rewritten.
It records the deck it read in deck.txt of the working directory, for the tests to compare.
"""
import os
import sys

name = os.path.basename(sys.argv[0]).lower()
role = os.environ.get("ENDF_STANDIN_ROLE", "clean")
here = os.getcwd()


def classic(x):
    if x == 0.0:
        return " 0.000000+0"
    m, e = f"{x:.6e}".split("e")
    s = f"{m}{'+' if int(e) >= 0 else '-'}{abs(int(e))}"
    return s.rjust(11)


def read_float(field):
    t = field.replace(" ", "")
    if not t:
        return 0.0
    if "e" in t.lower():
        return float(t)
    for i in range(len(t) - 1, 0, -1):
        if t[i] in "+-":
            return float(t[:i] + "e" + t[i:])
    return float(t)


def lossy(src, dst):
    out = []
    for line in open(src).read().splitlines():
        t = line.ljust(80)
        if t[70:72] == " 3" and t[72:75] != "  0":
            fields = [t[11 * k:11 * k + 11] for k in range(6)]
            new = []
            for k, f in enumerate(fields):
                try:
                    new.append(classic(read_float(f)) if f.strip() and ("." in f or k < 2) else f)
                except ValueError:
                    new.append(f)
            t = "".join(new) + t[66:]
        out.append(t[:75] + f"{len(out) + 1:5d}")
    open(dst, "w").write("\n".join(out) + "\n")


if name in ("fixup", "dictin", "endf2c"):
    deck = open(f"{name.upper()}.INP").read() if os.path.exists(f"{name.upper()}.INP") else ""
    open("deck.txt", "w").write(deck)
    if role == "crash":
        sys.stderr.write("Fortran runtime error: End of file\n")
        sys.exit(2)
    if role == "lossy":
        lossy("ENDFB.IN", "ENDFB.OUT")
    else:
        open("ENDFB.OUT", "w").write(open("ENDFB.IN").read())
    open(f"{name.upper()}.LST", "w").write(f" {name.upper()} stand-in listing\n")
    sys.exit(0)

lines = sys.stdin.read().splitlines()
open("deck.txt", "w").write("\n".join(lines) + "\n")
infile, outfile = lines[0], lines[1]
if not os.path.exists(infile):
    sys.exit(3)
if name == "stanef":
    if role == "crash":
        sys.stderr.write("Fortran runtime error: End of file\n")
        sys.exit(2)
    if role == "lossy":
        lossy(infile, outfile)
    else:
        open(outfile, "w").write(open(infile).read())
    print(" PROGRAM STANEF VERSION 8.10")
    sys.exit(0)

banner = {"checkr": "PROGRAM CHECKR VERSION 8.24", "fizcon": "PROGRAM FIZCON VERSION 8.19",
          "psyche": "PROGRAM PSYCHE VERSION 8.11", "inter": "PROGRAM INTER VERSION 8.08"}[name]
head = f"\n\n\n{banner}                              Run on 14-Nov-2026\n\nInput File Specification" \
       f"------------------------{infile}\nCheck the Entire File\n\n"
if role == "crash":
    open(outfile, "w").write(head)
    sys.stderr.write("At line 806 of file checkr.f (unit = 5, file = 'stdin')\nFortran runtime error: End of file\n")
    sys.exit(2)
if name == "psyche":
    body = ("TAPE BEING PROCESSED IS NUMBERED    0\n\n   CHECK MATERIAL 9999\n"
            "                   (NO WARNINGS DETECTED IN SECTIONS WITHOUT COMMENTS)\n" + "-" * 80 + "\nFILE  2\n"
            "   SECTION 151\n" + ("       AP NOT IN RANGE  1.11118E-01 TO  2.69858E-01\n" if role == "findings" else "")
            + "\n     CHECK ON ISOTOPE PROPERTIES\n\n" + "-" * 80 + "\n Done PSYCHE\n")
elif role == "findings":
    body = ("\nFILE IS UNSEQUENCED. NO SEQUENCE TESTS WILL BE DONE\n\n\n\n\n Check material   6-C - 12  MAT 9999\n"
            "\n\n  ERROR(S) FOUND IN MAT=9999, MF= 1, MT=451\n"
            "     LRP=2 ONLY ALLOWED IN DERIVED FILES              RECORD NUMBER       5\n"
            "\n\n  ERROR(S) FOUND IN MAT=9999, MF= 6, MT=  2\n"
            "     NE =   400 OUT OF RANGE     1 -   200            RECORD NUMBER     189\n"
            "     SECTION CANNOT BE CHECKED FROM SEQUENCE NUMBER      189 TO    2593\n"
            "\n     Encountered     3 errors,       0 warnings\n Done " + name.upper() + "\n")
else:
    body = "\nFILE IS UNSEQUENCED. NO SEQUENCE TESTS WILL BE DONE\n\n No problems to report\n\n Done " + \
           name.upper() + "\n"
open(outfile, "w").write(head + body)
print(f" {banner}\n")
sys.exit(0)
