#!/usr/bin/env python3
"""The command line that does what a script's answers to the old text menu did.

AZURE2 asked its run's questions in a text menu until the program took commands (``AZURE2 <command> [options]
FILE``); scripts answered them through a pipe, ``printf '3\\nn\\nparam.sav\\n\\n7\\n' | AZURE2 p.azr``
(their --no-readline is refused now: drop it).  This tool reads such answers, walks them through the menu's
questions as the menu read them (its order, its blank answers, its re-asking), and prints the command
and options that give the same run::

    python3 tools/menu_to_command.py '3\\nn\\nparam.sav\\n\\n7\\n'
    extrapolate --parameter-file=param.sav

The answers are given as ``printf`` would read them (``\\n`` between answers).  Options already on the old command
line that answered a question (``--covariance-band``) are given with ``--flags``.  Whether the menu asked for the
external-capture amplitude file depended on the project (a level with an external-capture multipolarity); the tool
tries both and says so when the two readings differ (``--capture`` / ``--no-capture`` choose one).  A file name is
taken as readable unless it is a number or a yes/no word (the menu asked again after a file it could not read).

Exit status: 0 with the command on standard output (notes on standard error); 1 when the answers do not make a run
(the menu's Exit, an end of input before a required answer, a y/n answer the menu refused): the reason on standard
error; 2 for a usage error.  The menu's table is the program's (src/RunRequest.cpp); this tool repeats its order
for the conversion of old scripts only and is checked against the program by tests/cli_commands/check.sh.
"""
import argparse
import re
import sys

# The menu's numbers and the commands (src/RunRequest.cpp, RunCommands()); 7 is Exit in a build with the sampler.
COMMANDS = {1: "calculate", 2: "fit", 3: "extrapolate", 4: "minos", 5: "rate", 6: "mcmc"}
EXIT_NUMBER = 7
YES = {"y", "yes", "1", "t", "true"}
NO = {"n", "no", "0", "f", "false"}


class NoRun(Exception):
    """The answers end the menu without a run."""


class Answers:
    """The answer lines as the menu reads them, one per question (getline)."""

    def __init__(self, lines):
        self.lines = list(lines)
        self.used = 0

    def next(self):
        if self.used >= len(self.lines):
            return None
        self.used += 1
        return self.lines[self.used - 1]


def leading_int(text):
    """operator>> of an int: the leading whole number of the text, or None."""
    m = re.match(r"\s*([+-]?\d+)", text)
    return int(m.group(1)) if m else None


def leading_float(text):
    m = re.match(r"\s*([+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?)", text)
    return float(m.group(1)) if m else None


def readable(name):
    """Whether the menu would have read `name` as a file (see the module's text)."""
    word = name.strip().lower()
    return bool(word) and leading_float(word) is None and word not in YES | NO


def simulate(lines, capture, band_given=False):
    """The command and options of the menu run the answers make; NoRun when they make none."""
    a = Answers(lines)
    notes = []
    while True:  # the menu's choice
        answer = a.next()
        if answer is None:
            raise NoRun("the input ends at the menu's choice (exit 1)")
        if not answer:
            continue
        n = leading_int(answer)
        if n is not None and 1 <= n <= EXIT_NUMBER:
            break
    if n == EXIT_NUMBER:
        raise NoRun("the answers choose Exit (no run, exit 0)")
    command = COMMANDS[n]
    options = []

    def flag(prompt):
        answer = a.next()
        if answer is None or not answer.strip():
            return False
        word = answer.strip().lower()
        if word in YES:
            return True
        if word in NO:
            return False
        raise NoRun(f'the menu refuses "{answer}" for "{prompt}" and stops (exit 1): the scripted-answer trap')

    def number(name, whole, low, high, low_open=False, not_zero=False, default=None):
        while True:
            answer = a.next()
            if answer is None:
                if default is None:
                    raise NoRun(f"the input ends before {name} (exit 1)")
                return default
            if default is not None and not answer.strip():
                return default
            x = leading_int(answer) if whole else leading_float(answer)
            if x is None:
                continue
            if (not_zero and x == 0) or (x <= low if low_open else x < low) or x > high:
                continue
            return x

    def word(name, required):
        while True:
            answer = a.next()
            if answer is None:
                if required:
                    raise NoRun(f"the input ends before {name} (exit 1)")
                return False
            if not required and not answer.strip():
                return False
            w = "".join(c for c in answer if c not in " \t\n")
            if w in ("yes", "no"):
                return w == "yes"

    def file(name, required):
        while True:
            answer = a.next()
            if answer is None:
                if required:
                    raise NoRun(f"the input ends before {name} (exit 1)")
                return ""
            answer = answer.rstrip(" \t")
            if not answer and not required:
                return ""
            if answer and readable(answer):
                return answer
            if answer:
                notes.append(f'"{answer}" is not taken as a file: the menu asked again')

    if command == "minos":
        options.append("--delta-chi2=%s" % fmt(number("the chi-squared variance", False, 0.0, float("inf"))))
    band = band_given
    if command in ("fit", "extrapolate", "minos") and not band_given:
        band = flag("Calculate cross-section uncertainty band?")
        if band:
            options.append("--covariance-band")
        if band and command in ("fit", "minos") and flag("Scale covariance to reduced chi-squared = 1?"):
            options.append("--scale-covariance")
    parameter = file("the parameter file", False)
    if parameter:
        options.append("--parameter-file=" + parameter)
    if command == "rate":
        while True:
            entrance = number("the rate's entrance pair", True, float("-inf"), float("inf"), not_zero=True)
            exit_pair = number("the rate's exit pair", True, -1, float("inf"), not_zero=True)
            if entrance != exit_pair:
                break
        options += ["--rate-entrance=%d" % entrance, "--rate-exit=%d" % exit_pair]
        if word("whether the temperatures come from a file", True):
            options.append("--temperature-file=" + file("the temperature file", True))
        else:
            for name in ("min", "max", "step"):
                options.append("--temperature-%s=%s" % (name, fmt(number("a temperature", False, 0.0, float("inf")))))
    if capture and command != "rate":
        amplitudes = file("the external-capture amplitude file", False)
        if amplitudes:
            options.append("--capture-amplitudes=" + amplitudes)
    if command == "mcmc":
        options.append("--walkers=%d" % number("the walkers", True, 2, 1000))
        options.append("--steps=%d" % number("the steps", True, 100, 10000000))
        spread = number("the parameter spread", False, 0.0, 50.0, low_open=True, default=5.0)
        energy = number("the energy spread", False, 0.0, 1.0, low_open=True, default=1.0)
        if spread != 5.0:
            options.append("--parameter-spread=" + fmt(spread))
        if energy != 1.0:
            options.append("--energy-spread=" + fmt(energy))
        options.append("--sampler-threads=%d" % number("the sampler's threads", True, 1, float("inf")))
        options.append("--basis=" + ("rwa" if word("the amplitude basis", True) else "physical"))
        if word("whether to overwrite samples.mcmc", False):
            options.append("--overwrite-samples")
    rest = a.lines[a.used:]
    if [line for line in rest if line.strip() not in ("", str(EXIT_NUMBER))]:
        notes.append("answers left unread: %r" % rest)
    return [command] + options, notes


def fmt(x):
    return ("%.12g" % x)


def lines_of(text):
    """The answer lines of a printf string: \\n escapes are line ends; a last line without its end still counts."""
    text = text.replace("\\n", "\n")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def translate(text, band_given=False, capture=None):
    """(the command line as a list, notes) for the answers `text`; NoRun when they make no run.  `capture`: whether
    the project has external capture (None: both readings, which must agree)."""
    lines = lines_of(text)
    readings = {}
    for c in ([capture] if capture is not None else [False, True]):
        try:
            readings[c] = simulate(lines, c, band_given)
        except NoRun as e:
            readings[c] = e
    runs = {c: r for c, r in readings.items() if not isinstance(r, NoRun)}
    if not runs:
        raise next(iter(readings.values()))
    if len({tuple(r[0]) for r in runs.values()}) > 1:
        raise NoRun("the answers mean different runs with and without external capture: %s; give --capture or "
                    "--no-capture" % " / ".join(" ".join(r[0]) for r in runs.values()))
    command, notes = next(iter(runs.values()))
    if len(runs) < len(readings):
        notes.append("valid only %s external capture" % ("with" if True in runs else "without"))
    return command, notes


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("answers", help="the menu answers as printf reads them, e.g. '2\\nn\\nparam.sav\\n\\n7\\n'")
    p.add_argument("--flags", default="", help="options of the old command line, e.g. '--covariance-band'")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--capture", dest="capture", action="store_const", const=True, default=None,
                   help="the project has external capture")
    g.add_argument("--no-capture", dest="capture", action="store_const", const=False)
    args = p.parse_args(argv)
    try:
        command, notes = translate(args.answers, "--covariance-band" in args.flags.split(), args.capture)
    except NoRun as e:
        print("menu_to_command: %s" % e, file=sys.stderr)
        return 1
    for n in notes:
        print("note: " + n, file=sys.stderr)
    print(" ".join(command))
    return 0


if __name__ == "__main__":
    sys.exit(main())
