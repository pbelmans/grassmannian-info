"""Preserve the existing presentation while migrating to static hosting.

Jinja runs only during export. Hugo owns the document shell and routes. Keeping
the content templates shared avoids dropping diagrams, prose, or interactions
while individual templates are migrated to Hugo in later commits.
"""

from datetime import datetime, timezone
from types import SimpleNamespace

from jinja2 import ChoiceLoader, DictLoader

from scripts.application import app, Dynkin, Grassmannian, grassmannians, horospherical_varieties
from scripts import literature


def environment():
    return app.jinja_env.overlay(loader=ChoiceLoader([
        DictLoader({"page.html": "{% block content %}{% endblock %}"}),
        app.jinja_loader,
    ]))


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
    env = environment()
    now = datetime.now(timezone.utc)
    return {route: env.get_template(template).render(now=now, **context)
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
