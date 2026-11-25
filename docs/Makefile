SPHINXBUILD ?= sphinx-build
SOURCEDIR = .
BUILDDIR = _build

.PHONY: help html pdf all

help:
	@$(SPHINXBUILD) --help

html:
	@$(SPHINXBUILD) -W --keep-going -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

pdf:
	@$(SPHINXBUILD) -W --keep-going -b pdf "$(SOURCEDIR)" "$(BUILDDIR)/pdf"

all: html pdf
