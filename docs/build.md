# Build notes

Python exports the mathematical data to `data/catalog.json`; Hugo renders the
pages, SVG diagrams and the root-system fallback from `layouts/`. The exporter
uses only the Python standard library. Sage and Oscar are needed only to
regenerate research data.

The catalogue and Hugo output are generated locally and in CI, and ignored by Git.
Exact invariants are exported as decimal strings; small ranks and node indices
remain JSON numbers.

## Verification

Build the site before running the tests. `scripts/verify.py` checks all 1,146
pages, exported invariants, citation targets, assets, aliases, Vogel point links,
certificate downloads and browser root-system calculations through rank 100.
The tests check the mathematical classifications and their presentation in the
built pages. Set `GRASSMANNIAN_PUBLIC_DIR` when testing a different output directory.

The workflow exports the data, builds Hugo, runs the verifier and tests, then deploys.
For a project-path build:

```shell
hugo --baseURL https://pbelmans.github.io/grassmannian-info/
python3 scripts/verify.py --base-url https://pbelmans.github.io/grassmannian-info/
```

## Root-system fallback

For classical root-system URLs beyond the prebuilt catalogue, such as `/A22`,
`404.html` uses the original root-description templates and browser arithmetic
to render the page, including its diagram, Cartan matrix and rank navigation.
Unknown routes retain the ordinary not-found message.

GitHub Pages serves this page with HTTP 404 even when the browser renders a
valid root system. The fallback covers root systems; Grassmannians outside the
stored data are not computed in the browser.
