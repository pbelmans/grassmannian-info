# Build notes

Python computes the invariants and renders the Jinja content in
`scripts/templates/` into `data/catalog.json`. Hugo supplies the document shell,
routes, assets and aliases. Flask is a build dependency; no Python server runs
on GitHub Pages. Sage and Oscar are only needed to regenerate research data.

The catalogue and Hugo output are generated locally and in CI, and ignored by Git.
Exact integers are exported as decimal strings to preserve their precision.

## Verification

`scripts/verify.py` compares all 1,146 generated page bodies with Flask
responses, including diagrams, controls, scripts and citations. It also checks
model data, aliases, Vogel point links, certificate downloads, and browser
root-system calculations through rank 100.

The workflow runs the unit tests, export, Hugo build and verifier before deployment.
For a project-path build:

```shell
hugo --baseURL https://pbelmans.github.io/grassmannian-info/
.venv/bin/python scripts/verify.py --base-url https://pbelmans.github.io/grassmannian-info/
```

## Root-system fallback

For classical root-system URLs beyond the prebuilt catalogue, such as `/A22`,
`404.html` uses the original root-description templates and browser arithmetic
to render the page, including its diagram, Cartan matrix and rank navigation.
Unknown routes retain the ordinary not-found message.

GitHub Pages serves this page with HTTP 404 even when the browser renders a
valid root system. The fallback covers root systems; Grassmannians outside the
stored data are not computed in the browser.
