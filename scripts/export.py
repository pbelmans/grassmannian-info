"""Export the existing mathematical models without porting their arithmetic.

Hugo renders the exported data; Python performs the mathematical computations.
"""

import json
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.geometry import Dynkin, Grassmannian, data, grassmannians, horospherical_varieties
from scripts import literature


def exact_json(value):
    """Keep exact invariants as strings; small node indices and plots stay numeric."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, dict):
        return {key: item if key in ("n", "k", "y", "z") and isinstance(item, int) else exact_json(item)
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [exact_json(item) for item in value]
    return value


def dynkin_record(D):
    return {
        "route": f"/{D.D}", "kind": "section", "family": "Dynkin type",
        "T": D.T, "n": D.n, "next_rank": D.n + 1, "previous_rank": D.n - 1, "type": D.D, "lie": D.lie(),
        "invariants": {
            "dimension": D.dimension(), "rank": D.n,
            "Weyl group order": D.weyl_group_cardinality(),
            "Coxeter number": D.coxeter_number(),
            "exponents": D.exponents(), "roots": D.number_of_roots(),
        },
        "cartan": D.cartan_matrix(),
        "positive_roots": D.number_of_positive_roots(),
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
        "name": G.plaintext(), "heading": G.name(), "html": G.html(), "latex": G.latex(),
        "type": G.D.D, "k": G.k,
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
        "latex": X.latex(), "heading": X.name(), "html": X.html(),
        "type": X.D.D, "y": X.y, "z": X.z, "orbits": [X.Y.url(), X.Z.url()],
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
    target = ROOT / "data"
    target.mkdir(exist_ok=True)
    type_names = {row["type"] for row in exported if row["family"] == "Dynkin type"}
    type_names.update(row["automorphisms"]["type"] for row in exported if row["family"] == "Grassmannian")
    types = {D: dynkin_record(Dynkin(D[0], int(D[1:]))) for D in sorted(type_names)}
    # Low-rank presentations used by Freudenthal and the table's isomorphisms.
    extra = {}
    for T, n, k in (("B", 1, 1), ("C", 1, 1), ("D", 1, 1), ("D", 2, 1), ("D", 2, 2)):
        G = Grassmannian(T, n, k)
        row = {"isomorphisms": G.isomorphisms(), "route": "/" + G.url()}
        if (T, n) != ("D", 1):
            row.update(latex=G.latex(), invariants={"dimension": G.dimension(), "index": G.index(), "Euler characteristic": G.rank()})
        extra[row["route"]] = row
    payload = {"records": {row["route"]: row for row in exported},
               "types": exact_json(types), "extra": exact_json(extra),
               "ranks": {T: {str(n): sorted(parabolics) for n, parabolics in ranks.items()}
                         for T, ranks in grassmannians.items()},
               "horospherical_types": {f"{X.D.T}-{X.D.n}-{X.y}-{X.z}": "/horospherical/" + X.plaintext() for X in horospherical_varieties.values()},
               "reference_order": [ref["key"] for ref in literature.ordered_references(literature.REFERENCES)]}
    path = target / "catalog.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n")
    counts = {family: sum(row["family"] == family for row in exported)
              for family in sorted({row["family"] for row in exported})}
    print(json.dumps({"seconds": round(perf_counter() - started, 3),
                      "bytes": path.stat().st_size, "pages": counts}))


if __name__ == "__main__":
    main()
