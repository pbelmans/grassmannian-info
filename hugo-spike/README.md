# Hugo migration experiment

The existing site can be served statically. Python computes the catalogue at
build time; Sage is only needed to regenerate the stored research data. There
is no user-specific server state.

Feature preservation is a requirement. The initial minimal Hugo templates
omitted substantial parts of the site and have been replaced. This branch now
uses the existing Jinja content templates during the export, and Hugo supplies
the document shell, routes, assets, aliases, and output. **This is a hybrid build,
not a completed conversion to native Hugo templates.** No Flask server runs in
production for the prebuilt pages. Vogel work is outside this experiment.

## Build and verify

From the repository root, with the Python dependencies and Hugo installed:

```sh
.venv/bin/python hugo-spike/export.py
hugo --source hugo-spike
.venv/bin/python hugo-spike/verify.py
hugo server --source hugo-spike --port 1317
```

In a separate worktree, use the main checkout's Python executable if needed.
Node is required for the verifier's browser-arithmetic checks. Regenerate the
export after editing Python, Jinja, or mathematical source data: Hugo's watcher
does not run Python. Generated JSON and HTML are ignored.

One measured complete export took 15.2 seconds and the Hugo build 344 ms with
Hugo 0.161.1. The generated catalogue, including page content, is about 32 MB.

## Presentation coverage

| Pages or features | Preserved in the static preview |
| --- | --- |
| Home page | Original periodic table, classifications, colours, information toggle, hover panel, and cross-highlighting |
| 862 Grassmannians | Marked clickable Dynkin diagrams, alternative realisations, vertical Betti sidebar, automorphisms, projective geometry and duality, exceptional constructions, quantum spectra and multiplicities, visual Lefschetz grids, HPD, higher Fano properties, and citations |
| 191 horospherical varieties | Original overview and colour controls, orbit data, marked diagrams, blowups and projections, exceptional collections, quantum spectra, and HPD boxes |
| 81 prebuilt root systems | Original root descriptions, Weyl data, Cartan matrices, diagrams, and rank navigation |
| Freudenthal | Both squares, marked diagrams, Lie algebras, all five tabs, and the information toggle |
| Bibliography | Overview, chronological entries, topic filter, published links, arXiv and MR identifiers, and summaries |
| Other pages and assets | About, explanation, site navigation, website icons, MathJax/Xy-pic, and quantum-certificate download |
| Legacy aliases | B1/1 and C1/1 redirect to A1/1 |

`verify.py` compares all 1,140 generated page bodies with independent Flask
responses, including markup, text, SVG, scripts, and citation targets. This
would reject the earlier abbreviated layouts. It also round-trips all model
records, checks the aliases and certificate bytes, and compares browser
root-system calculations with Python in 88 cases through rank 100. Exact
integers remain strings across JSON/Hugo: 10,427 invariant entries exceed
`2**53 - 1`, and the largest has 255 digits.

Safari checks cover the C6/3 diagram, vertical Betti numbers, spectrum and HPD
grid; the horospherical colours and information toggle; and Freudenthal's
expanded values. The Freudenthal placeholder zeroes and two incorrect family
formulas were corrected in a separate commit, also applied to Flask.

## Unbounded root-system URLs

Flask generates classical root-system pages beyond the stored catalogue, such
as `/A21`. A finite export cannot enumerate all such URLs. The static fallback
uses the original root-description templates and a small browser calculation
with `BigInt` for Weyl orders. It retains diagrams, Cartan matrices and rank
navigation. Ordinary catalogue pages do not need this calculation.

Static hosting must serve `404.html` for otherwise unmatched URLs for this
fallback to run, or rewrite classical type URLs to that document with status
200. The preview uses Hugo's 404 handling. Production routing must be checked
before switching hosts.

The old site also links to Grassmannians outside its stored projective data
(for example A21/1 and D21/21); these already fail in Flask. They are separate
from the working root-system URLs above. The comparison is a regression check,
not a claim that all historical mathematical data and links are correct.

## Finishing the native Hugo conversion

Convert content templates one page family at a time, keeping the complete
Flask comparison as a gate. The current compatibility layer deliberately keeps
the original presentation instead of maintaining a second incomplete design.
Extract the mathematical models from `application.py` to remove Flask as a
build dependency. Keep exact polynomial arithmetic in Python, and retain the
stored Sage results. Only then remove the corresponding Jinja content templates.

This resembles the hyperkaehler.info rewrite's data preparation step more than
Fanography's use of Go template helpers for enrichment. The content adapter
provides the Hugo routes; the substantial arithmetic should not be ported to Go
templates. Production deployment and CI have not been switched.
