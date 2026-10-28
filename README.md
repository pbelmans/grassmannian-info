# Grassmannian.info

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21296785.svg)](https://doi.org/10.5281/zenodo.21296785)

A periodic table of (generalised) Grassmannians: [grassmannian.info](https://grassmannian.info).

A static Hugo site. Python computes the mathematical data at build time; Hugo
renders the pages. The legacy `flask` branch is no longer maintained.

## Structure

- `hugo.toml` — Hugo configuration.
- `content/` — the content adapter generates the pages from the catalogue.
- `layouts/` — Hugo page templates, diagrams and shared components.
- `data/` — Sage source data and the generated catalogue.
- `static/` — CSS, JavaScript and images, published under `/static/`.
- `scripts/` — mathematical models in `geometry.py`, data export and verification.
- `references/` — bibliography and stored mathematical results.

## Building

Install Hugo 0.161.1 and Python 3.13 or later, then:

```shell
python3 scripts/export.py
hugo server --port 1317
```

Regenerate the export after changing Python or source data. To build and
verify the production output (Node is required for the browser checks):

```shell
hugo
python3 -m unittest discover -s tests -v
python3 scripts/verify.py
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
