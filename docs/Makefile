# Minimal Sphinx makefile.
#   make -C docs html      -> docs/_build/html/index.html
#   make -C docs api       -> docs/api/html/index.html (Doxygen, C++ internals)
SPHINXOPTS  ?=
SPHINXBUILD ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = _build

.PHONY: help html api clean
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

api:
	@doxygen Doxyfile

clean:
	@rm -rf "$(BUILDDIR)" api
