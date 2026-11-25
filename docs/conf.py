"""Sphinx configuration for the airnet user and API documentation."""

from __future__ import annotations

import sys
from pathlib import Path

from docutils import nodes

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from airnet.__about__ import __version__

project = "airnet"
copyright = "2023-present, Oak Ridge National Laboratory and contributors"
version = release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "rst2pdf.pdfbuilder",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_preserve_defaults = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
primary_domain = "py"
nitpicky = True
nitpick_ignore = [
    ("py:class", "collections.abc.Mapping"),
    ("py:class", "enum.Enum"),
    ("py:class", "optional"),
]

html_theme = "pydata_sphinx_theme"
html_title = f"airnet {release} documentation"
html_theme_options = {
    "show_toc_level": 2,
    "navigation_with_keys": True,
}

pdf_documents = [
    # rst2pdf requires a four-item tuple even when no author is credited.
    ("index", "airnet", f"airnet {release} Documentation", ""),
]
pdf_stylesheets = ["sphinx", "airnet"]
pdf_use_index = True
pdf_use_modindex = True
pdf_fit_mode = "shrink"


def _replace_pdf_abbreviations(app, doctree, docname):
    """Downgrade Sphinx signature abbreviations unsupported by rst2pdf."""
    if app.builder.name != "pdf":
        return
    for abbreviation in list(doctree.findall(nodes.abbreviation)):
        abbreviation.replace_self(nodes.inline("", abbreviation.astext()))


def setup(app):
    """Register documentation-build compatibility hooks."""
    app.connect("doctree-resolved", _replace_pdf_abbreviations)
