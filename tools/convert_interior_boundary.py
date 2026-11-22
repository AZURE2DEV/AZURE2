#!/usr/bin/env python3
"""Rewrite a project written for the branch feature/upgrade-6-interior-potential (potential records with
`interior_boundary`, `l` and `interior_charge_e1/e2`) in the `wave` layout of docs/upgrade/potential_background.md.

    python3 tools/convert_interior_boundary.py old.azr new.azr [--lmax N] [--final-l L[,L...]]

Only the <potentials> section changes.  The records of one pair become one record:

  old                                              new
  ---------------------------------------------    ------------------------------------------------------------
  a record without `l`                             the pair's record (its central term, radii, step)
  a record with `l 1 3`, depth D                   `wave l 1 depth D MeV ...`, `wave l 3 depth D MeV ...`
  `interior_boundary b`                            `... background interior boundary b capture interior`
  `interior_charge_e1 q` (`_e2`)                   `... charge_e1 q` (`charge_e2 q`) on those waves
  a pair with records for single l only            the first is the pair's record; every l keeps its depth as a wave

What cannot be carried over is refused with the reason, and nothing is written:

  - records of one pair that differ in more than the depth (a Woods-Saxon radius or diffuseness, the Coulomb radius,
    the matching radius, the step, a spin-orbit term): a `wave` has its own depth (and Gaussian width) only;
  - a depth that is not positive (the `wave` layout has attractive wells only).

Two things change in the physics, and the conversion says so on stderr:

  - a partial wave that holds an external-capture FINAL state gets no background (a final state in a background wave
    is refused in the `wave` layout): it keeps its depth, and its channel functions are those of the potential
    outside the radius alone.  The final-state waves are read from <levels> (a capture row's final J^pi and energy,
    then that level's particle rows); `--final-l` names them instead.
  - a bound final state must be a bound state of its wave's potential within 1 keV of the level's energy, and the
    interior is normalised to the level's ANC; the engine refuses the project otherwise and prints the depth that
    binds there.

A record with `interior_boundary` and no `l` acted on every partial wave: it becomes one `wave` per l from 0 to
`maxOrbitalMomentum` of <config>, or to `--lmax`.
"""
import argparse
import re
import sys


def fail(message):
    sys.exit("convert_interior_boundary: " + message)


def section(lines, name):
    """(first, last) line indices of the section's body, or None."""
    try:
        a = next(i for i, line in enumerate(lines) if line.strip() == f"<{name}>")
        b = next(i for i, line in enumerate(lines) if line.strip() == f"</{name}>" and i > a)
    except StopIteration:
        return None
    return a + 1, b


def final_waves(lines):
    """The orbital angular momenta of the particle channels of every external-capture final state."""
    body = section(lines, "levels")
    rows = [line.split() for line in lines[body[0]:body[1]]] if body else []
    rows = [r for r in rows if len(r) >= 27]
    finals = {(r[14], r[15], float(r[16])) for r in rows if r[26] == "10" and r[12] != "0"}
    waves = set()
    for r in rows:
        if r[26] == "0" and any(r[0] == j and r[1] == pi and abs(float(r[2]) - e) < 1e-3 for j, pi, e in finals):
            waves.add(int(r[7]) // 2)
    return waves


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("old")
    parser.add_argument("new")
    parser.add_argument("--lmax", type=int)
    parser.add_argument("--final-l", help="comma-separated l of the final-state partial waves (default: from <levels>)")
    args = parser.parse_args()
    lines = open(args.old).read().split("\n")
    body = section(lines, "potentials")
    if not body:
        fail("no <potentials> section")
    records, current = [], None
    for line in lines[body[0]:body[1]]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if not line[0].isspace():
            words = line.split()
            if len(words) != 2 or words[0] != "potential":
                fail(f"not a potential record: '{line}'")
            current = {"name": words[1], "fields": []}
            records.append(current)
        elif current is None:
            fail(f"a field before the first record: '{line}'")
        else:
            key, _, value = line.strip().partition(" ")
            current["fields"].append((key, value.strip()))
    old_words = ("interior_boundary", "interior_charge_e1", "interior_charge_e2", "l")
    if not any(k in old_words for r in records for k, _ in r["fields"]):
        fail("no record has interior_boundary, interior_charge_e1/e2 or l: nothing to convert")
    if any(k == "wave" for r in records for k, _ in r["fields"]):
        fail("a record already has a wave field")
    config = section(lines, "config")
    lmax = args.lmax
    if lmax is None and config:
        for line in lines[config[0]:config[1]]:
            m = re.match(r"\s*maxOrbitalMomentum\s*=\s*(\d+)", line)
            if m:
                lmax = int(m.group(1))
    finals = {int(v) for v in args.final_l.split(",")} if args.final_l else final_waves(lines)

    def get(record, key):
        values = [v for k, v in record["fields"] if k == key]
        if len(values) > 1:
            fail(f"potential {record['name']}: {key} is given twice")
        return values[0] if values else None

    def central(record):
        words = (get(record, "central") or "").split()
        if not words or words[0] not in ("woods_saxon", "gaussian"):
            fail(f"potential {record['name']}: no central woods_saxon or gaussian term")
        parts = {words[i]: (words[i + 1], words[i + 2]) for i in range(1, len(words) - 2, 3)}
        if float(parts["depth"][0]) <= 0.0:
            fail(f"potential {record['name']}: depth {parts['depth'][0]} is not positive; the wave layout has attractive wells only")
        return words[0], parts

    pairs = {}
    for r in records:
        pairs.setdefault((get(r, "pair"), get(r, "pair_key")), []).append(r)
    out, notes = [], []
    for (pair, _), group in pairs.items():
        whole = [r for r in group if get(r, "l") is None]
        if len(whole) > 1:
            fail(f"pair {pair}: more than one record without l")
        base = whole[0] if whole else group[0]
        base_shape, base_parts = central(base)
        waves = []
        for r in group:
            shape, parts = central(r)
            for key in ("coulomb_radius", "matching_radius", "step", "spin_orbit"):
                if get(r, key) != get(base, key):
                    fail(f"pair {pair}: the records {base['name']} and {r['name']} differ in {key} "
                         f"('{get(base, key)}' and '{get(r, key)}'); a wave has its own depth only")
            if shape != base_shape:
                fail(f"pair {pair}: the records {base['name']} and {r['name']} have different central shapes")
            for key in ("radius", "diffuseness"):
                if shape == "woods_saxon" and parts[key] != base_parts[key]:
                    fail(f"pair {pair}: the records {base['name']} and {r['name']} differ in the Woods-Saxon {key} "
                         f"({' '.join(base_parts[key])} and {' '.join(parts[key])}); a wave has its own depth only")
            boundary = get(r, "interior_boundary")
            listed = get(r, "l")
            if listed is None:
                if boundary is None:
                    continue
                if lmax is None:
                    fail(f"potential {r['name']}: interior_boundary without l acts on every partial wave; "
                         "give --lmax (no maxOrbitalMomentum in <config>)")
                ls = list(range(lmax + 1))
            else:
                ls = [int(v) for v in listed.split()]
            for l in ls:
                text = f"l {l}"
                if r is not base:
                    text += f" depth {' '.join(parts['depth'])}"
                    if shape == "gaussian" and parts["width"] != base_parts["width"]:
                        text += f" width {' '.join(parts['width'])}"
                if boundary is not None and l in finals:
                    notes.append(f"pair {pair}, l = {l}: holds an external-capture final state, so it gets no background "
                                 f"(it had interior_boundary {boundary}); a bound final state must be the potential's own within 1 keV")
                elif boundary is not None:
                    text += f" background interior boundary {boundary} capture interior"
                    for word, new in (("interior_charge_e1", "charge_e1"), ("interior_charge_e2", "charge_e2")):
                        if get(r, word) is not None and float(get(r, word)) != 0.0:
                            text += f" {new} {get(r, word)}"
                if text != f"l {l}":
                    waves.append((l, text))
        seen = [l for l, _ in waves]
        if len(seen) != len(set(seen)):
            fail(f"pair {pair}: an l is named by more than one record")
        out.append(f"potential {base['name']}")
        for key, value in base["fields"]:
            if key in old_words:
                continue
            out.append(f"    {key:<15} {value}")
            if key == "central":
                out.extend(f"    {'wave':<15} {text}" for _, text in sorted(waves))
        if len(group) > 1:
            notes.append(f"pair {pair}: the records {', '.join(r['name'] for r in group)} are now the one record {base['name']}")
    new = lines[:body[0]] + out + lines[body[1]:]
    with open(args.new, "w") as f:
        f.write("\n".join(new))
    for note in notes:
        print("convert_interior_boundary: " + note, file=sys.stderr)
    print(f"convert_interior_boundary: wrote {args.new} ({sum(1 for line in out if ' wave ' in line)} wave fields)", file=sys.stderr)


if __name__ == "__main__":
    main()
