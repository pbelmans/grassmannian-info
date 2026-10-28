# Grassmannian.info

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21296785.svg)](https://doi.org/10.5281/zenodo.21296785)

A periodic table of (generalised) Grassmannians: [grassmannian.info](https://grassmannian.info).

A static Hugo site, with mathematical computations and Jinja rendering in Python
at build time. The legacy `flask` branch is no longer maintained.

## Structure

- `hugo.toml` — Hugo configuration.
- `content/` — the content adapter generates the pages from the catalogue.
- `layouts/` — Hugo templates.
- `data/` — Sage source data and the generated catalogue.
- `static/` — CSS, JavaScript and images, published under `/static/`.
- `scripts/` — mathematical models in `geometry.py`, export, verification and Jinja content templates.
- `references/` — bibliography and stored mathematical results.

## Building

Install Hugo 0.161.1 and Python, then:

```shell
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/export.py
hugo server --port 1317
```

Regenerate the export after changing Python, Jinja or source data. To build and
verify the production output (Node is required for the browser checks):

```shell
.venv/bin/python -m unittest discover -s tests -v
hugo
.venv/bin/python scripts/verify.py
```

After changing Sass, rebuild the stylesheet:

```shell
sass --no-source-map assets/style.scss static/style.css
```

See [research data](docs/research.md) for Sage and Oscar computations, and
[build notes](docs/build.md) for verification and the 404 fallback.

## Deployment

Pushes to `main` build, verify and deploy `public/` to GitHub Pages through
`.github/workflows/pages.yml`.

## How to cite

The website is archived on Zenodo. To cite the version-independent latest
release, use the concept DOI
[10.5281/zenodo.21296785](https://doi.org/10.5281/zenodo.21296785).

```bibtex
@online{grassmannian,
  author = {Belmans, Pieter},
  title  = {Grassmannian.info --- A periodic table of (generalised) Grassmannians},
  url    = {https://grassmannian.info},
  doi    = {10.5281/zenodo.21296785},
  year   = {2026},
}
```
