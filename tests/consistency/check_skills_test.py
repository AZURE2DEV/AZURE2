#!/usr/bin/env python3
"""Unit test of tools/check_skills.py: each rule of docs/upgrade/RULES.md §3 it enforces
fails on a minimal skill tree that breaks it, and a clean tree passes.

    python3 tests/consistency/check_skills_test.py
"""

import datetime as dt
import os
import sys
import tempfile
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import check_skills  # noqa: E402

TODAY = dt.date(2026, 11, 6)
GOOD_DESC = "Do one thing with AZURE2. Use for X and Y; for Z use other-skill instead."


class SkillTree:
    def __init__(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        os.makedirs(os.path.join(self.root, "docs", "source"))
        with open(os.path.join(self.root, "docs", "source", "page.rst"), "w") as f:
            f.write("Page\n====\n")

    def write(self, rel, text):
        path = os.path.join(self.root, ".claude", "skills", rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        return path

    def skill(self, name, body="Body.\n", desc=GOOD_DESC, front_name=None):
        return self.write(f"{name}/SKILL.md", f"---\nname: {front_name or name}\ndescription: {desc}\n---\n\n{body}")

    def run(self, user_dir=None, strict=False):
        return check_skills.Checker(self.root, TODAY, strict=strict, user_skills_dir=user_dir).run()

    def close(self):
        self.tmp.cleanup()


class CheckSkillsTest(unittest.TestCase):
    def setUp(self):
        self.t = SkillTree()

    def tearDown(self):
        self.t.close()

    def rules(self, findings, errors_only=True):
        return sorted({f.rule for f in findings if not (errors_only and f.warning)})

    def test_clean_tree_passes(self):
        self.t.skill("azure2-a", "See `azure2-b` and docs/source/page.rst and [ref](references/x.md).\n"
                                 "Builds before 2026-10-27 (expires: 2027-02-01) wrote 6 digits.\n")
        self.t.write("azure2-a/references/x.md", "detail\n")
        self.t.skill("azure2-b", desc="Another task entirely. Use for W; for X use azure2-a.")
        self.assertEqual(self.t.run(), [])

    def test_r4_sizes(self):
        self.t.skill("s", "x\n" * 600)
        f = self.t.run()
        self.assertEqual(self.rules(f), [])           # 500 < n <= 800 is a warning
        self.assertEqual(self.rules(f, errors_only=False), ["R4"])
        self.assertEqual(self.rules(self.t.run(strict=True)), ["R4"])
        self.t.skill("s", "x\n" * 801)
        self.assertEqual(self.rules(self.t.run()), ["R4"])
        self.t.skill("s")
        self.t.write("s/references/long.md", "y\n" * 301)
        self.assertEqual(self.rules(self.t.run()), ["R4"])

    def test_r5_history(self):
        for line in ("Fixed by 3fba1192 in the engine.",
                     "The run (job 1420848) died.",
                     "Behaviour since dev changed.",
                     "We saw this on 2026-10-05.",
                     "Builds before 2026-10-27 differ.",                        # no expiry
                     "Builds before 2026-10-27 (expires: 2026-11-01) differ.",  # expired
                     "Builds before 2026-10-27 (expires: 2027-06-01) differ."): # beyond 90 days
            with self.subTest(line=line):
                self.t.skill("s", line + "\n")
                self.assertEqual(self.rules(self.t.run()), ["R5"])

    def test_r5_allow_comment_and_non_hashes(self):
        self.t.skill("s", "Example of a forbidden note: fixed in dev 3fba1192 <!-- check_skills: allow R5 -->\n"
                          "Numbers like 2.6e-10, deadbeef words and 1234567 are not hashes.\n")
        self.assertEqual(self.t.run(), [])

    def test_r6_names_and_owners(self):
        self.t.skill("s", front_name="other")
        self.assertEqual(self.rules(self.t.run()), ["R6"])
        self.t.skill("s")
        self.t.skill("t")                                   # same description as s
        self.assertEqual(self.rules(self.t.run()), ["R6"])
        self.t.skill("t", "Use `azure2-gone` for that.\n", desc="Unique text.")
        self.assertEqual(self.rules(self.t.run()), ["R6"])

    def test_r6_user_level_collision(self):
        self.t.skill("s")
        with tempfile.TemporaryDirectory() as user:
            os.makedirs(os.path.join(user, "s"))
            with open(os.path.join(user, "s", "SKILL.md"), "w") as f:
                f.write("---\nname: s\ndescription: something else\n---\n")
            self.assertEqual(self.rules(self.t.run(user_dir=user)), ["R6"])

    def test_r7_description_length(self):
        self.t.skill("s", desc="x" * 601)
        self.assertEqual(self.rules(self.t.run()), ["R7"])

    def test_r8_personal_paths(self):
        for line in ("Run /groups/someone/R-matrix/rmfit.", "Open /Users/me/AZURE2.",
                     "cd ~/AZURE2/build", "See /home/me/x."):
            with self.subTest(line=line):
                self.t.skill("s", line + "\n")
                self.assertEqual(self.rules(self.t.run()), ["R8"])
        self.t.skill("s", "The table is in ~/.azure2/ and $RMFIT_HOME/rmfit; awk 's!~/^ *0/' and t[i]~/x/.\n")
        self.assertEqual(self.t.run(), [])

    def test_r8_scripts(self):
        self.t.skill("s")
        self.t.write("s/scripts/run.py", "DATA = '/groups/x/y'\n")
        self.assertEqual(self.rules(self.t.run()), ["R8"])

    def test_links(self):
        for line in ("See docs/source/missing.rst.", "See references/none.md.", "[x](../nowhere/SKILL.md)"):
            with self.subTest(line=line):
                self.t.skill("s", line + "\n")
                self.assertEqual(self.rules(self.t.run()), ["links"])

    def test_links_to_repository_files(self):
        os.makedirs(os.path.join(self.t.root, "scripts"))
        open(os.path.join(self.t.root, "scripts", "tool.py"), "w").close()
        self.t.skill("s", "Run scripts/tool.py.\n")
        self.assertEqual(self.t.run(), [])

    def test_main_exit_status(self):
        self.t.skill("s")
        self.assertEqual(check_skills.main(["--root", self.t.root, "--today", "2026-11-06", "--user-skills", "none"]), 0)
        self.t.skill("s", "x\n" * 900)
        self.assertEqual(check_skills.main(["--root", self.t.root, "--today", "2026-11-06", "--user-skills", "none"]), 1)


if __name__ == "__main__":
    unittest.main()
