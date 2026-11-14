"""A test helper (not part of pyazr): writes small ENDF-6 files by hand for the tests of pyazr.endf.

Each record is written field by field as ENDF-102 (2023) lays it out (pp. 30-36); the directory of MF=1/MT=451 and
the sequence numbers are made from the sections (pp. 33, 41-42).  The numbers are written in the form asked
(``pyazr.endf.NUMBER_FORMS``: classic, e, widest), so that a test can read back what it wrote.
"""

from pyazr.endf import NUMBER_FORMS


class Writer:
    """Records of one section as lines of 66 columns of content (MAT, MF, MT, NS are added by :func:`material`)."""

    def __init__(self, form="classic"):
        self.f = NUMBER_FORMS[form]
        self.rows = []

    def _i(self, n):
        return f"{n:11d}"

    def cont(self, c1, c2, l1, l2, n1, n2):
        self.rows.append(self.f(c1) + self.f(c2) + "".join(self._i(v) for v in (l1, l2, n1, n2)))
        return self

    def text(self, t):
        self.rows.append(t.ljust(66)[:66])
        return self

    def list(self, c1, c2, l1, l2, n2, values):
        self.cont(c1, c2, l1, l2, len(values), n2)
        for k in range(0, len(values), 6):
            self.rows.append("".join(self.f(v) for v in values[k:k + 6]))
        return self

    def _ranges(self, nbt, laws):
        pairs = [v for ab in zip(nbt, laws) for v in ab]
        for k in range(0, len(pairs), 6):
            self.rows.append("".join(self._i(v) for v in pairs[k:k + 6]))

    def tab1(self, c1, c2, l1, l2, nbt, laws, x, y):
        self.cont(c1, c2, l1, l2, len(nbt), len(x))
        self._ranges(nbt, laws)
        xy = [v for ab in zip(x, y) for v in ab]
        for k in range(0, len(xy), 6):
            self.rows.append("".join(self.f(v) for v in xy[k:k + 6]))
        return self

    def tab2(self, c1, c2, l1, l2, nbt, laws, nz):
        self.cont(c1, c2, l1, l2, len(nbt), nz)
        self._ranges(nbt, laws)
        return self


def material(path, sections, mat=9999, za=6012.0, awr=11.896910, awi=0.998623, nsub=10010, lrp=-1, nmod=1,
             text=("  6-C - 12 TEST       EVAL-NOV26 pyazr tests", " hand-written           DIST-NOV26",
                   "---- TEST                MATERIAL 9999", "----- INCIDENT PROTON DATA", "------ ENDF-6"),
             form="classic", sequence=True, tpid=" a file of the pyazr tests", emax=2.0e7):
    """Write a material: ``sections`` {(MF, MT): Writer} in any order; MF=1/MT=451 is made here (header, text,
    directory with NC the lines of each section without SEND, p. 41)."""
    f = NUMBER_FORMS[form]
    keys = sorted(sections)
    nxc = len(keys) + 1
    head = Writer(form)
    head.cont(za, awr, lrp, 0, 0, nmod).cont(0.0, 0.0, 0, 0, 0, 6).cont(awi, emax, 0, 0, nsub, 8)
    head.cont(0.0, 0.0, 0, 0, len(text), nxc)
    for t in text:
        head.text(t)
    n451 = 4 + len(text) + nxc
    for mf, mt in [(1, 451)] + keys:
        nc = n451 if (mf, mt) == (1, 451) else len(sections[(mf, mt)].rows)
        head.rows.append(" " * 22 + "".join(f"{v:11d}" for v in (mf, mt, nc, 0)))
    allsec = dict(sections)
    allsec[(1, 451)] = head
    out = [tpid.ljust(66) + f"{1:4d}{0:2d}{0:3d}{0:5d}" if sequence else tpid.ljust(66) + f"{1:4d}{0:2d}{0:3d}"]
    zero = f(0.0) * 2 + f"{0:11d}" * 4

    def line(content, mf, mt, ns):
        return content.ljust(66) + f"{mat:4d}{mf:2d}{mt:3d}" + (f"{ns:5d}" if sequence else "")

    last_mf = None
    for mf, mt in sorted(allsec):
        if last_mf is not None and mf != last_mf:
            out.append(zero + f"{mat:4d}{0:2d}{0:3d}" + (f"{0:5d}" if sequence else ""))
        for k, row in enumerate(allsec[(mf, mt)].rows, start=1):
            out.append(line(row, mf, mt, k))
        out.append(line(zero, mf, 0, 99999))
        last_mf = mf
    out.append(zero + f"{mat:4d}{0:2d}{0:3d}" + (f"{0:5d}" if sequence else ""))
    out.append(zero + f"{0:4d}{0:2d}{0:3d}" + (f"{0:5d}" if sequence else ""))
    out.append(zero + f"{-1:4d}{0:2d}{0:3d}" + (f"{0:5d}" if sequence else ""))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return path
