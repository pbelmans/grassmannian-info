# Hugo feasibility prototype

A static rewrite is feasible. Keep the mathematical models in Python, evaluate
them at build time, and let Hugo render their results. The expensive Sage
computations already produce stored data; ordinary page requests do not run
Sage. There is no user-specific server state that needs to survive the rewrite.

This directory is an isolated experiment, not the completed rewrite. It leaves
the Flask application, deployment workflow, and Vogel work unchanged.

## Reproduce

Use the repository's Python environment (`pip install -r requirements.txt`)
and Hugo. Tested with Hugo 0.161.1. From the repository root:

```sh
.venv/bin/python hugo-spike/export.py
hugo --source hugo-spike
.venv/bin/python hugo-spike/verify.py
hugo server --source hugo-spike --port 1317
```

In a separate worktree, the Python executable can come from the main checkout's
virtual environment. The exporter imports the models from its own worktree.
Generated JSON and HTML are ignored. Regenerate the data after changing a
mathematical model or source dataset; `hugo server` alone does not run Python.

## Measured result

On the development machine, one complete export took 4.73 seconds and a Hugo
build took 147 milliseconds. These are measurements of the minimal prototype
templates, not a promise about a finished site's build time.

| Exported records | Count |
| --- | ---: |
| Grassmannians | 862 |
| Horospherical varieties | 191 |
| Dynkin types | 81 |
| Total model pages checked | 1,134 |

Hugo writes 1,137 pages including the home page, horospherical index, and
bibliography. The generated JSON is 6.21 MB; the output occupies about 12 MB.
The verifier independently obtains records from the existing Python models and
checks every record after its passage through JSON, Hugo, and an HTML parser.
It also checks the displayed invariants and local citation targets on every
model page. This verifies the export boundary, not the underlying mathematical
theorems or visual equivalence with the current site.

## Where the computations belong

| Work | Static architecture |
| --- | --- |
| Sage generation of degrees, Hilbert values, and spectra | Keep the stored results and separate regeneration scripts |
| Exact Betti polynomials, Weyl orders, dimensions, indices | Evaluate the existing Python models once per build |
| Quantum classifications, Lefschetz shapes, HPD and citations | Evaluate `literature.py` at build time |
| Page structure, bibliography, Dynkin SVG diagrams | Hugo templates and partials |
| Filtering, highlighting, eigenvalue plots | Existing browser JavaScript, adapted to generated data |

The exporter temporarily imports `application.py`, so Flask remains a build
dependency of this experiment. A production rewrite should extract `Dynkin`,
`Grassmannian`, `Horospherical`, and the polynomial helpers into a module without
Flask. The same models can serve both implementations during comparison. There
is no reason to port exact polynomial arithmetic into Go templates.

## Comparison with the earlier rewrites

The local Fanography rewrite computes derived properties in
`layouts/partials/lib/enrich-fano.html`. That is convenient for its data, but
Grassmannian.info has substantial exact polynomial arithmetic worth retaining.

The local hyperkaehler.info rewrite is the closer precedent:
`scripts/build-data.mjs` prepares data before Hugo runs and serializes large
integers as strings; research outputs remain separate. Follow that architecture
here, retaining Python as the computation language. Both sites use content
adapters to generate pages from data, which this prototype also demonstrates.
See the [Hugo content adapter documentation](https://gohugo.io/content-management/content-adapters/).

## Remaining migration work

1. Extract the models and define a documented export schema. Preserve exact
   integers as strings: 10,427 scalar invariant entries in this export exceed
   `2**53 - 1`, and the longest has 255 digits. Plot coordinates can remain
   floating-point numbers. Keep structured data in `hugo.Data`; the prototype
   passes only a lookup key through page parameters, avoiding Hugo's
   normalization of parameter-map keys.
2. Port the Jinja templates and macros to Hugo. There are about 2,400 template
   lines excluding Vogel, including the tables, Dynkin diagrams, and
   Freudenthal page. Preserve the existing MathJax/Xy-pic setup initially;
   changing the math renderer is a separate migration.
3. Define the route manifest. The prototype preserves case and generates
   `/C6/3/` and `/horospherical/X3(3,2)/`. Static hosting must serve the existing
   slashless links through its directory redirects. Flask also accepts type
   pages beyond the finite dataset, such as `/A21`, and some isomorphism links
   point beyond that dataset. Choose and enforce a supported range, or generate
   additional type pages. Include legacy aliases and the quantum-certificate
   download. This prototype does not implement those aliases, downloads, prose
   pages, diagrams, or interactive controls.
4. Compare both renderings across the route manifest, including citations,
   unknown versus negative results, diagrams, spectra, and mobile layouts.
   The prototype stores more model data than it displays; a full port must
   render all existing boxes and status information.
5. Switch CI to generate data and run Hugo, then deploy static files. Neither
   Python nor Sage would be needed on the serving host. Update navigation to
   match the supported route manifest before switching production.

The computation is a good fit for a short build step. Template conversion and
checking existing URLs and interaction account for most of the remaining work.
