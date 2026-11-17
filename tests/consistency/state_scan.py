#!/usr/bin/env python3
"""Find the state of the engine that outlives a call: thread_local variables, non-const static and namespace-scope
variables, mutable members (docs/upgrade/evaluation.md, the state table; check 38 of repo_consistency_test.py).

A source-level scan, not a compiler: comments, string and character literals are blanked, braces are matched to
tell a namespace, a class and a block apart, and each declaration is classified by its words.  What it reports:

  thread_local  any variable declared thread_local (every scope)
  static        a variable declared static in a block (function-level) or a static data member of a class,
                when it is not const (a pointer to const is not const: `static const T *p` is reported)
  namespace     a variable at namespace scope (static or not, in a named or an unnamed namespace) that is not const
  mutable       a member declared mutable
  once          a std::once_flag member of a class: the member it guards is written once, in a const method, behind
                a pointer or not (CoulombTable::ContentIdentity), which the kinds above do not see

std::once_flag, std::mutex and std::atomic objects are reported like any other variable: a mutex guards state, a
flag or counter is state.  Run as a script it prints one line per finding: file, symbol, kind, line.

    python3 tests/consistency/state_scan.py [directory ...]
"""
import os
import re
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
# The engine's sources (not external/, not the tests).
SCANNED = ("src", "include", "api", "physics", "experiment", "stopping", "convert", "io", "support", "inference")
SKIPPED = ("/tests/", "/build", "/_deps/")
SUFFIXES = (".cpp", ".cc", ".h", ".hpp")


def blank(text):
    """The text with comments, string and character literals replaced by blanks (line breaks kept)."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(re.sub(r"[^\n]", " ", text[i:j]))
            i = j
        elif c == "R" and text.startswith('R"', i) and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == "_")):
            m = re.match(r'R"([^(\s]*)\(', text[i:])
            if m:
                end = text.find(")" + m.group(1) + '"', i)
                j = n if end < 0 else end + len(m.group(1)) + 2
                out.append('"' + re.sub(r"[^\n]", " ", text[i + 1:j - 1]) + '"')
                i = j
            else:
                out.append(c)
                i += 1
        elif c in "\"'":
            if c == "'" and i > 0 and (text[i - 1].isalnum()):  # digit separator 1'000
                out.append(c)
                i += 1
                continue
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            out.append(c + re.sub(r"[^\n]", " ", text[i + 1:j - 1]) + c)
            i = j
        else:
            out.append(c)
            i += 1
    text = "".join(out)
    # preprocessor lines (with their continuations)
    return re.sub(r"(?m)^[ \t]*#(?:[^\n]*\\\n)*[^\n]*", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


KEYWORD_HEAD = re.compile(r"^(using|typedef|template|friend|extern|class|struct|union|enum|namespace|static_assert|"
                          r"return|if|else|for|while|do|switch|case|default|goto|break|continue|throw|delete|new|"
                          r"public|private|protected|try|catch|co_return|operator)\b")


def scope_of(head):
    """The kind of scope a '{' opens, from the text before it."""
    h = head.strip()
    if re.search(r"\bnamespace\b[\w\s:]*$", h) or re.search(r'\bextern\s+"\s*"\s*$', h):
        return "namespace"
    if re.search(r"\b(class|struct|union)\b[^(){};=]*$", h) and not re.search(r"\b(class|struct|union)\s+\w+\s+\w+\s*$", h):
        return "class"
    if re.search(r"\benum\b[^(){};]*$", h):
        return "enum"
    if h.endswith("=") or h.endswith(",") or h.endswith("(") or h.endswith("[") or h.endswith("{") or \
            h.endswith("return") or re.search(r"[\w>\]]\s*$", h) and re.search(r"[=]\s*[\w:<>]*\s*$", h):
        return "init"
    if re.search(r"\)\s*(const|noexcept|override|final|mutable|->[^{]*)*\s*$", h) or h.endswith("try") or \
            h.endswith("do") or h.endswith("else") or h == "" or h.endswith(";") or h.endswith(":") or h.endswith("}"):
        return "block"
    # a declarator with a brace initialiser: T name{...}
    if re.search(r"[\w>\]]\s*$", h):
        return "init"
    return "block"


def declaration(stmt, scope="namespace"):
    """(name, words) if the statement declares variables, else None.  `stmt` has no initialiser braces."""
    s = " ".join(stmt.split())
    if not s or KEYWORD_HEAD.match(s):
        return None
    head = re.split(r"(?<![=!<>])=(?!=)", s, maxsplit=1)[0]
    head = re.sub(r"\{.*", "", head).strip()
    if "(" in head:
        # a function declaration, a call, or a variable with a parenthesised initialiser (T x(args)): only the
        # last counts, recognised by a known object type in front
        m = re.match(r"^(?:(?:static|thread_local|inline|mutable|volatile)\s+)*(std::atomic\b[^()]*?)\s+(\w+)\s*\(", head)
        if not m and scope == "block":
            # in a block, `static T name(args)` and `thread_local T name(args)` can only be variables
            m = re.match(r"^((?:(?:static|thread_local|inline|volatile)\s+)+[^()]*?)\s*\b(\w+)\s*\(", head)
            if m and not re.search(r"\b(static|thread_local)\b", m.group(1)):
                m = None
        if not m:
            return None
        # more declarators with their own parenthesised initialisers: T a(0), b(0)
        more = re.findall(r"\)\s*,\s*(\w+)\s*\(", head[m.end(2):])
        return ",".join([m.group(2)] + more), head[:m.start(2)]
    # several declarators (T a, *b, c[2]): split at the commas outside template brackets; each takes the type words of
    # the first, plus its own '*' or '&'
    pieces, depth, cur = [], 0, ""
    for ch in head:
        depth += (ch == "<") - (ch == ">")
        if ch == "," and depth == 0:
            pieces.append(cur)
            cur = ""
        else:
            cur += ch
    pieces.append(cur)
    m = re.search(r"([A-Za-z_]\w*)\s*(\[[^\]]*\]\s*)*$", pieces[0])
    if not m:
        return None
    name, words = m.group(1), pieces[0][:m.start(1)]
    others = []
    for p in pieces[1:]:
        q = re.match(r"^\s*([*&\s]*(?:const\s*)?)([A-Za-z_]\w*)\s*(\[[^\]]*\]\s*)*$", p)
        if not q:
            return None
        others.append(q.group(2))
    if others:
        name = ",".join([name] + others)
    if not re.search(r"[\w>*&]\s*$", words) or re.match(r"^\s*$", words):
        return None
    if re.match(r"^(?:const|constexpr|static|inline|volatile|mutable|thread_local|unsigned|signed|long|short)\s*$",
                words.strip()) and words.strip() in ("static", "inline", "mutable", "thread_local"):
        return None
    return name, words


def is_const(words):
    """True when the declared object itself is const (a pointer to const is not)."""
    w = words.strip()
    while re.search(r"<[^<>]*>", w):  # template arguments do not qualify the object
        w = re.sub(r"<[^<>]*>", "", w)
    if re.search(r"\bconstexpr\b", w):
        return True
    if "*" in w or "&" in w:
        last = re.split(r"[*&]", w)[-1]
        if "&" in w and "*" not in w:  # a reference cannot be reseated; what it refers to decides
            return bool(re.search(r"\bconst\b", w))
        return bool(re.search(r"\bconst\b", last))
    return bool(re.search(r"\bconst\b", w))


def scan_text(text):
    """[(symbol, kind, line)] for one source text."""
    text = blank(text)
    found = []
    stack = []  # scope kinds
    start = 0
    depth_init = 0
    i, n = 0, len(text)
    # the statements are cut at ';', '{' and '}', tracking the scope
    while i < n:
        c = text[i]
        if c == "{":
            kind = scope_of(text[start:i]) if depth_init == 0 else "init"
            if kind == "init":
                depth_init += 1
                stack.append("init")
            else:
                stack.append(kind)
                start = i + 1
        elif c == "}":
            kind = stack.pop() if stack else "block"
            if kind == "init":
                depth_init -= 1
            else:
                start = i + 1
        elif c == ";" and depth_init == 0:
            stmt = text[start:i]
            scope = next((k for k in reversed(stack) if k != "init"), "namespace")
            line = text.count("\n", 0, start + len(stmt) - len(stmt.lstrip())) + 1
            found += classify(stmt, scope, line)
            start = i + 1
        elif c == "(":
            # skip a parenthesised group whole (a for-header's ';', a lambda's braces inside an argument list are
            # handled by the brace logic when they appear outside parentheses only)
            d, j = 1, i + 1
            while j < n and d:
                if text[j] == "(":
                    d += 1
                elif text[j] == ")":
                    d -= 1
                elif text[j] in "{}" and d:
                    break
                j += 1
            if d == 0:
                i = j
                continue
        i += 1
    return found


def classify(stmt, scope, line):
    s = " ".join(stmt.split())
    # a label or an access specifier in front ("public:", "case 1:")
    s = re.sub(r"^(?:(?:public|private|protected)\s*:\s*)+", "", s)
    if scope == "enum" or not s:
        return []
    tl = bool(re.search(r"\bthread_local\b", s))
    st = bool(re.search(r"\bstatic\b", s))
    mu = bool(re.search(r"(?<!\])\s*\bmutable\b", s)) and scope == "class"
    once = scope == "class" and not st and not mu and bool(re.match(r"^std::once_flag\s+\w+\s*$", s))
    if once:
        return [(s.split()[-1], "once", line)]
    if not (tl or st or mu or scope == "namespace"):
        return []
    d = declaration(s, scope)
    if not d:
        return []
    name, words = d
    if re.search(r"\b(class|struct|enum|using|typedef)\b", words):
        return []
    names = name.split(",")
    if tl:
        return [(n, "thread_local", line) for n in names]
    if mu:
        return [(n, "mutable", line) for n in names]
    if is_const(words):
        return []
    if scope == "namespace":
        if re.search(r"\bextern\b", words):
            return []
        return [(n, "namespace", line) for n in names]
    if st and scope in ("block", "class"):
        return [(n, "static", line) for n in names]
    return []


def scan(root=ROOT, directories=SCANNED):
    """{(file, symbol, kind): [lines]} for the scanned tree."""
    out = {}
    for d in directories:
        for dirpath, _, files in os.walk(os.path.join(root, d)):
            for f in sorted(files):
                if not f.endswith(SUFFIXES):
                    continue
                path = os.path.join(dirpath, f)
                rel = os.path.relpath(path, root)
                if any(x in "/" + rel for x in SKIPPED):
                    continue
                with open(path, encoding="utf-8", errors="replace") as fh:
                    for symbol, kind, line in scan_text(fh.read()):
                        out.setdefault((rel, symbol, kind), []).append(line)
    return out


def read_allowed(path):
    """{(file, symbol, kind): reason} from allowed_state.tsv (tab-separated; '#' starts a comment line)."""
    allowed = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4 or not cols[3].strip():
                allowed[("?", raw.strip(), "?")] = ""
                continue
            allowed[(cols[0], cols[1], cols[2])] = cols[3]
    return allowed


if __name__ == "__main__":
    found = scan(directories=tuple(sys.argv[1:]) or SCANNED)
    for (rel, symbol, kind), lines in sorted(found.items()):
        print(f"{rel}\t{symbol}\t{kind}\t{','.join(map(str, lines))}")
    counts = {}
    for (_, _, kind) in found:
        counts[kind] = counts.get(kind, 0) + 1
    print("# " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), file=sys.stderr)
