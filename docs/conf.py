"""Sphinx configuration for the opticalRadiation + radiationDose docs."""

project = "opticalRadiation + radiationDose"
author = "DeGroot Research Group"
copyright = "2026, DeGroot Research Group"

extensions = [
    "myst_parser",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "colon_fence",
]

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

html_theme = "furo"
html_title = "opticalRadiation + radiationDose"
html_static_path = ["_static"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

mathjax3_config = {
    "tex": {
        "macros": {
            "RR": "{\\mathbb{R}}",
            "vec": ["{\\mathbf{#1}}", 1],
        }
    }
}
