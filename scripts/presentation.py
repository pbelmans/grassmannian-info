"""Render page content with Jinja; Hugo supplies the document shell and routes."""

from datetime import datetime, timezone
from functools import cache
from pathlib import Path
from types import SimpleNamespace

from jinja2 import Environment, FileSystemLoader, select_autoescape

from scripts.geometry import Dynkin, Grassmannian, Horospherical, grassmannians, horospherical_varieties
from scripts import literature


@cache
def environment():
    env = Environment(
        loader=FileSystemLoader(Path(__file__).parent / "templates"),
        autoescape=select_autoescape(),
    )
    env.globals.update(
        Dynkin=Dynkin, Grassmannian=Grassmannian, Horospherical=Horospherical,
        reference_catalog=literature.REFERENCES,
        reference_aliases={
            alias: key for key, ref in literature.REFERENCES.items()
            for alias in ref.get("aliases", [])
        },
    )
    return env


def render_template(template, **context):
    return environment().get_template(template).render(now=datetime.now(timezone.utc), **context)


def page_specs(model_records):
    for record in model_records:
        route = record["route"]
        if record["family"] == "Dynkin type":
            D = route[1:]
            yield route, "type.html", {"D": Dynkin(D[0], int(D[1:]))}
        elif record["family"] == "Grassmannian":
            D, k = route[1:].split("/")
            yield route, "grassmannian.html", {"G": Grassmannian(D[0], int(D[1:]), int(k))}
        else:
            yield route, "horospherical.html", {"X": horospherical_varieties[record["name"]]}
    yield "/", "table.html", {"grassmannians": grassmannians}
    yield "/horospherical", "horospherical.table.html", {"horospherical": horospherical_varieties}
    yield "/freudenthal", "freudenthal.html", {}
    yield "/vogel", "vogel.html", {}
    yield "/about", "about.html", {}
    yield "/explained", "explained.html", {}
    yield "/projective-duality", "projective-duality.html", {"grassmannians": grassmannians,
        "references": literature.ordered_references(["tevelev", "magic-square", "kuznetsov-hyperplanes"])}
    yield "/references", "references.html", {"references": literature.ordered_references(literature.REFERENCES)}


def render_pages(model_records):
    return {route: render_template(template, **context)
            for route, template, context in page_specs(model_records)}


def type_templates():
    """Reuse the original formulas for classical types beyond the build range."""
    class Rank:
        def __str__(self):
            return "__RANK__"

        def __index__(self):
            # The browser supplies the Cartan entries, without a rank bound.
            return 0

        def __add__(self, other):
            assert other == 1
            return "__RANK_PLUS_1__"

        def __sub__(self, other):
            assert other == 1
            return "__RANK_MINUS_1__"

    env = environment()
    result = {}
    for T in "ABCD":
        module = env.get_template("type.html").make_module({"D": Dynkin(T, 4)})
        D = SimpleNamespace(
            T=T, n=Rank(), dimension=lambda: "__DIMENSION__",
            coxeter_number=lambda: "__COXETER__", exponents=lambda: ["__EXPONENTS__"],
            number_of_roots=lambda: "__ROOTS__",
            weyl_group_cardinality=lambda: "__WEYL_ORDER__",
        )
        cartan = str(module.cartan(D)).replace("\\begin{pmatrix}", "\\begin{equation}\\begin{pmatrix}__CARTAN__")
        cartan = cartan.replace("\\end{pmatrix}", "\\end{pmatrix}\\end{equation}")
        result[T] = {"basic": str(module.basic(D)), "roots": str(module.roots(D)), "cartan": cartan}
    return result
