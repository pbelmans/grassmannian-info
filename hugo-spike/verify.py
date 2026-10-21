"""Compare complete Hugo page content with Flask, not just exported numbers."""

from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess

from export import exact_json, records
from application import app, Dynkin


class Page(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.model = ""
        self.ids, self.citations = set(), set()
        self.tokens = []
        self.in_model = False
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.ids.add(attributes["id"])
        if attributes.get("href", "").startswith("#reference-"):
            self.citations.add(attributes["href"][1:])
        if attributes.get("id") == "model-data":
            self.in_model = True
        self.tokens.append((tag, sorted(attrs)))

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        self.handle_endtag(tag)

    def handle_endtag(self, tag):
        if tag == "script":
            self.in_model = False
        self.tokens.append(("/" + tag,))

    def handle_data(self, value):
        if self.in_model:
            self.model += value
        value = " ".join(value.split())
        if value:
            self.tokens.append(("text", value))


def content(source):
    source = source.split('<div class="container" style="max-width: 1200px">', 1)[1]
    source = source.rsplit("</div>", 1)[0]
    return re.sub(r'<script type="application/json" id="model-data">.*?</script>', "", source, flags=re.S)


def compare(route, public, client):
    source = (public / route.lstrip("/") / "index.html").read_text()
    response = client.get(route)
    assert response.status_code == 200, (route, response.status_code)
    expected = Page(content(response.get_data(as_text=True))).tokens
    actual = Page(content(source)).tokens
    assert actual == expected, (route, "page content differs from Flask")
    page = Page(source)
    assert not page.citations - page.ids, (route, page.citations - page.ids)
    return page


def verify_dynamic_types():
    cases = [(T, n) for T, minimum in (("A", 1), ("B", 2), ("C", 2), ("D", 3))
             for n in (*range(minimum, 22), 37, 100)]
    script = """
      const fs = require('node:fs');
      const {model} = require('./hugo-spike/static/static/dynkin-page.js');
      const cases = JSON.parse(fs.readFileSync(0, 'utf8'));
      process.stdout.write(JSON.stringify(cases.map(([T, n]) => model(T, n))));
    """
    result = subprocess.run(["node", "-e", script], input=json.dumps(cases),
                            text=True, capture_output=True, check=True,
                            cwd=Path(__file__).parents[1])
    for (T, n), actual in zip(cases, json.loads(result.stdout), strict=True):
        D = Dynkin(T, n)
        assert actual == {
            "dimension": str(D.dimension()), "roots": str(D.number_of_roots()),
            "order": str(D.weyl_group_cardinality()), "coxeter": str(D.coxeter_number()),
            "exponents": D.exponents(), "cartan": D.cartan_matrix(),
        }, (T, n)
    print(f"Verified browser root-system arithmetic against Python in {len(cases)} cases, through rank 100.")


def main():
    public = Path(__file__).parent / "public"
    count, longest, large = 0, 0, 0
    client = app.test_client()
    for record in records():
        page = compare(record["route"], public, client)
        assert json.loads(page.model) == exact_json(record), record["route"]
        for value in record["invariants"].values():
            for number in value if isinstance(value, list) else [value]:
                if isinstance(number, int):
                    longest = max(longest, len(str(number)))
                    large += abs(number) > 2**53 - 1
        count += 1
    for route in ("/", "/horospherical", "/freudenthal", "/projective-duality", "/references", "/about", "/explained"):
        compare(route, public, client)
        count += 1
    for alias in ("B1/1", "C1/1"):
        source = (public / alias / "index.html").read_text()
        assert 'http-equiv="refresh"' in source and "/A1/1/" in source
    certificates = "references/quantum-certificates.json"
    assert (public / certificates).read_bytes() == (Path(__file__).parents[1] / certificates).read_bytes()
    print(f"Verified {count} complete page bodies against Flask, including diagrams, controls, scripts, and citations.")
    print(f"Checked model JSON, aliases, and the quantum-certificate download.")
    print(f"Longest integer: {longest} digits; {large} invariant values exceed 2**53 - 1.")
    verify_dynamic_types()


if __name__ == "__main__":
    main()
