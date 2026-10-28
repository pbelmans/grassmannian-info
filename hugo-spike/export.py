"""Export the existing mathematical models without porting their arithmetic.

This prototype imports application.py, so Flask is a build dependency for now.
A production migration should first move the models to a Flask-free module.
"""

import json
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from application import Dynkin, Grassmannian, data, horospherical_varieties
import literature
from presentation import render_pages, type_templates


def exact_json(value):
    """All exact integers travel as decimal strings; plot floats stay numeric."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, dict):
        return {key: exact_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [exact_json(item) for item in value]
    return value


def dynkin_record(D):
    return {
        "route": f"/{D.D}", "kind": "section", "family": "Dynkin type",
        "invariants": {
            "dimension": D.dimension(), "rank": D.n,
            "Weyl group order": D.weyl_group_cardinality(),
            "Coxeter number": D.coxeter_number(),
            "exponents": D.exponents(), "roots": D.number_of_roots(),
        },
        "cartan": D.cartan_matrix(),
    }


def common_record(X):
    points = X.eigenvalues()
    return {
        "kind": "page",
        "invariants": {
            "dimension": X.dimension(), "index": X.index(),
            "Euler characteristic": X.rank(), "Betti numbers": X.betti(),
        },
        "eigenvalues": list(points) if points is not None else [],
    }


def grassmannian_record(G):
    record = common_record(G)
    record.update({
        "route": "/" + G.url(), "family": "Grassmannian",
        "name": G.plaintext(), "latex": G.latex(),
        "levi": [D.D for D in G.L],
        "automorphisms": {
            "type": G.Aut().D, "dimension": G.Aut().dimension(),
            "components": G.pizero(),
        },
        "isomorphisms": sorted(G.isomorphisms()),
        "flags": {name: getattr(G, "is_" + name)()
                  for name in ("minuscule", "cominuscule", "adjoint", "coadjoint")},
        "small_quantum": G.small_quantum_result(),
        "big_quantum": G.big_qh_semisimple(),
        "exceptional": G.exceptional_sequences(),
        "lefschetz": G.lefschetz(), "hpd": G.hpd_results(),
        "projective_duality": G.projective_duality(),
        "references": G.references(), "fano": G.fano(),
    })
    record["invariants"].update({
        "embedding dimension": G.embedding(), "degree": G.degree(),
        "Hilbert values": G.hilbert_series(),
    })
    return record


def horospherical_record(X):
    record = common_record(X)
    # These are the existing horospherical classifications, without new claims.
    small = ("semisimple" if X.small_qh_semisimple() else
             "non-semisimple" if X.small_qh_not_semisimple() else "unknown")
    record.update({
        "route": "/horospherical/" + X.plaintext(),
        "family": "Horospherical variety", "name": X.plaintext(),
        "latex": X.latex(), "orbits": [X.Y.url(), X.Z.url()],
        "codimensions": [X.codimY(), X.codimZ()],
        "small_quantum": {"status": small},
        "big_quantum": X.big_qh_semisimple(),
        "exceptional": X.exceptional_sequences(), "lefschetz": X.lefschetz(),
        "fano": X.fano(),
    })
    return record


def records():
    # Vogel also links to low-rank types and A21 outside the Grassmannian data.
    types = sorted({row["type"] for row in data} | {"A21", "B1", "C1", "D2"},
                   key=lambda D: (D[0], int(D[1:])))
    for D in types:
        yield dynkin_record(Dynkin(D[0], int(D[1:])))
    for row in data:
        D = row["type"]
        yield grassmannian_record(Grassmannian(D[0], int(D[1:]), row["parabolic"]))
    for X in horospherical_varieties.values():
        yield horospherical_record(X)


def main():
    started = perf_counter()
    exported = exact_json(list(records()))
    target = Path(__file__).parent / "data"
    target.mkdir(exist_ok=True)
    payload = {"records": {row["route"]: row for row in exported},
               "pages": render_pages(exported),
               "type_templates": type_templates(),
               "references": exact_json(literature.ordered_references(literature.REFERENCES))}
    path = target / "catalog.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n")
    counts = {family: sum(row["family"] == family for row in exported)
              for family in sorted({row["family"] for row in exported})}
    print(json.dumps({"seconds": round(perf_counter() - started, 3),
                      "bytes": path.stat().st_size, "pages": counts}))


if __name__ == "__main__":
    main()
