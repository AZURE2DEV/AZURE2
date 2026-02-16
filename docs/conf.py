# Configuration file for the Sphinx documentation builder.

project = 'AZURE2'
copyright = '2015-2026, E. Uberseder and R. J. deBoer'
author = 'E. Uberseder and R. J. deBoer'
release = '2.0'

extensions = [
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_theme_options = {
    'logo_only': False,
    'version_selector': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False,
}

html_title = 'AZURE2 Documentation'
html_short_title = 'AZURE2'

# LaTeX output configuration
latex_elements = {
    'preamble': r'\usepackage{amsmath}',
}
