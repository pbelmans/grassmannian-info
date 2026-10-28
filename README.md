# Grassmannian.info

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21296785.svg)](https://doi.org/10.5281/zenodo.21296785)

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

## How to deploy the website

With GitHub CLI installed and configured, simply run

```shell
gh workflow run deploy.yml --ref flask
```

## Local development

```shell
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
sh script
```

The compiled stylesheet is included for deployment. After changing Sass, rebuild it with:

```shell
sass --no-source-map assets/style.scss static/style.css
```

## Research references and checks

`references/catalog.json` contains the bibliography used by `/references` and
the individual Grassmannian pages. Its primary links point to published
versions where available, with arXiv and MathSciNet identifiers retained.
Retrieve new bibliographical metadata with `bibgetter pull arXiv:...` or
`bibgetter pull MR...`, then update the catalog's publication details and
summary. The mathematical hypotheses and applicability rules live in
`literature.py`; unknown results must remain distinct from negative results.

Run the page, citation, and mathematical consistency checks with:

```shell
.venv/bin/python -m unittest discover -s tests -v
```

The additional small-quantum semisimplicity cases have exact characteristic
polynomials in `references/quantum-certificates.json`. Coefficients are in
ascending order, for multiplication by the ample generator at `q=1` over `QQ`.
SageMath can verify the stored square-free polynomials and independent known
examples, or regenerate the matrices and certificates:

```shell
sage scripts/quantum_semisimplicity.sage --check
sage scripts/quantum_semisimplicity.sage
```

Generated intermediates and local research drafts are ignored.

Projective-dual degrees through Dynkin rank 8 are stored in
`references/projective-duality.json`. The Julia generator uses **Oscar's
IntersectionTheory** (`tn_variety`, `tn_bundle`, `chern_class`, and `integral`)
to evaluate Tevelev's Theorems 7.1–7.2 (arXiv version). It computes the Taylor
expansion at `q = 1` of `q ∫ c_top(Ω_X ⊗ O_X(q))`; the first nonconstant
coefficient is the degree of the dual, and its exponent is the codimension.
Truncated power series keep the computation practical even for `E8/4`.

The generator uses exact rational arithmetic, verifies the degree of the
original Grassmannian, and repeats each computation with two regular
cocharacters. For defective varieties it also checks the vanishing of the
higher jet Chern numbers. The values use the ample generator, including for
spinor varieties. Empty duals have null degree and codimension.

The isolated Julia environment pins Oscar 1.8.2 and does not require Julia
at runtime, either for Flask or for the Hugo export.

```shell
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts scripts/projective_duality.jl
julia --project=scripts scripts/projective_duality.jl --check
julia --project=scripts scripts/projective_duality.jl --case C6/3
```

An independent Sage implementation checks the same numbers by directly
expanding the jet-bundle Chern polynomial, using integer arithmetic:

```shell
sage -python scripts/projective_duality.py --check
```
