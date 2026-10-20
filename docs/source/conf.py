"""Sphinx configuration for the AZURE2 documentation.

Reconstructed: the repository tracked only the generated HTML under
``docs/_build``, never the sources or this file, so the documentation could not
be rebuilt. The .rst files were recovered from ``_build/html/_sources`` and the
theme and extensions inferred from the generated output. Adjust freely -- this
is a starting point that reproduces the previous look, not a historical record.
"""

project = "AZURE2"
copyright = "University of Notre Dame"
author = "AZURE2 developers"

# Keep in step with AZURE2_VERSION in the top-level CMakeLists.txt.
# The version comes from the CMakeLists.txt version block, the single source.
import os, re
_cmake = open(os.path.join(os.path.dirname(__file__), "..", "..", "CMakeLists.txt")).read()
_v = {k: re.search(r"set\(AZURE2_VERSION_%s (\d+)\)" % k, _cmake).group(1) for k in ("MAJOR", "MINOR", "PATCH")}
release = "%s.%s.%s" % (_v["MAJOR"], _v["MINOR"], _v["PATCH"])
version = "%s.%s" % (_v["MAJOR"], _v["MINOR"])

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "AZURE2 Documentation"

# Do not copy the prompt when someone copies a shell example.
copybutton_prompt_text = r"\$ |>>> |\.\.\. "
copybutton_prompt_is_regexp = True
