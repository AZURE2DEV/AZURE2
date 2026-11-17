#!/usr/bin/env python3
"""A project written with partial widths -> the same project written with a gamma_decay record, by hand: the
definitions of docs/upgrade/gamma_parametrisation.md (Eqs. 2.1, 2.2, 2.4), not the engine's map.  The radiative widths
of the level in <levels> become 0 and the record gives Gamma_gamma = sum of them, b_f = Gamma_f / Gamma_gamma and,
for a final state with two multipolarities, delta_f = -sign(Gamma_L Gamma_L+1) (Gamma_L+1 / Gamma_L)^1/2
(Krane-Steffen); a final state whose widths are all negative gets 'sign -'.

    python3 make_record.py IN.azr OUT.azr LEVEL_ID REMAINDER_E_f [fit|fixed]

REMAINDER_E_f: the final state's E_x (as <levels> writes it) whose branching is the remainder.  Every number of the
record is written with 17 significant digits (repr), so the map's products are the only rounding.
"""
import math
import sys

src, dst, level, remainder = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
mode = sys.argv[5] if len(sys.argv) > 5 else "fit"
lines = open(src, encoding="utf-8", newline="").read().split("\n")
widths = {}  # E_f text -> {L: width}
jpi = None
out = []
inside = False
for l in lines:
    if l.startswith("<levels>"):
        inside = True
    elif l.startswith("</levels>"):
        inside = False
    t = l.split()
    if inside and len(t) > 30 and t[8] == level and t[26] == "10":
        widths.setdefault(t[16], {})[int(t[7]) // 2] = float(t[11])
        J = float(t[0])
        jpi = (f"{int(2 * J)}/2" if J != int(J) else str(int(J))) + ("+" if t[1] == "1" else "-")
        t[11] = "0"
        l = " " + " ".join(t)
    out.append(l)
total = sum(sum(abs(w) for w in ws.values()) for ws in widths.values())
rec = ["<gammaDecays>", f"gamma_decay level{level}", f"    level        {level}", f"    J            {jpi}",
       f"    total_width  {total!r} eV {mode}"]
mixing, signs = [], []
for ef, ws in widths.items():
    label = "g+15O" + ("" if float(ef) == 0.0 else f"[{ef}]")
    b = sum(abs(w) for w in ws.values()) / total
    if b == 0.0 and len(ws) == 1:
        continue  # a final state not listed: b = 0, fixed
    if ef == remainder:
        rec.append(f"    remainder    {label}")
    else:
        rec.append(f"    branching    {label} {b!r} {mode}")
    Ls = sorted(ws)
    first = ws[Ls[0]] if ws[Ls[0]] != 0.0 else (ws[Ls[1]] if len(Ls) > 1 else 0.0)
    if first < 0.0:
        signs.append(f"    sign         {label} -")
    if len(Ls) == 2 and ws[Ls[1]] != 0.0:
        gl, gl1 = ws[Ls[0]], ws[Ls[1]]
        d = -math.copysign(1.0, gl * gl1) * math.sqrt(abs(gl1) / abs(gl))
        mixing.append(f"    mixing       {label} {d!r} {mode}")
rec += mixing + signs + ["</gammaDecays>", ""]
text = "\n".join(out)
if not text.endswith("\n"):
    text += "\n"
open(dst, "w", encoding="utf-8", newline="\n").write(text + "\n".join(rec))  # LF on every platform
