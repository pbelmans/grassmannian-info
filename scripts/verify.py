"""Verify generated pages, mathematical data, assets and browser computations."""

from html.parser import HTMLParser
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.export import exact_json, records
from scripts.geometry import Dynkin
from scripts.presentation import render_pages


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


def check_site_links(source, site_root):
    for link in re.findall(r'(?:href|src)="(/[^/"][^"]*|/)"', source):
        assert link.startswith(site_root), (link, "link outside site path")
    for link in re.findall(r'\\href\{(/[^} ]+)', source):
        assert link.startswith(site_root), (link, "math link outside site path")
    if site_root != "/":
        assert '"xlink:href", "/"' not in source, "unprefixed Dynkin node links"


def compare(route, public, expected_html, site_root):
    source = (public / route.lstrip("/") / "index.html").read_text()
    check_site_links(source, site_root)
    expected = Page(expected_html).tokens
    # Strip only URL prefixes, not the repository name inside external links.
    normalized = re.sub(r'(["{])' + re.escape(site_root), r'\1/', content(source))
    actual = Page(normalized).tokens
    assert actual == expected, (route, "page content differs from Jinja rendering")
    page = Page(source)
    assert not page.citations - page.ids, (route, page.citations - page.ids)
    return page


def verify_dynamic_types():
    cases = [(T, n) for T, minimum in (("A", 1), ("B", 2), ("C", 2), ("D", 3))
             for n in (*range(minimum, 22), 37, 100)]
    script = """
      const fs = require('node:fs');
      const {model} = require('./static/dynkin-page.js');
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


def verify_fallback(public, site_root):
    script = r'''
      const assert = require('node:assert/strict');
      const fs = require('node:fs');
      const vm = require('node:vm');
      const {publicDir, root} = JSON.parse(fs.readFileSync(0, 'utf8'));
      const html = fs.readFileSync(publicDir + '/404.html', 'utf8');
      const config = html.match(/<script type="text\/x-mathjax-config">([\s\S]*?)<\/script>/)[1];
      let settings;
      vm.runInNewContext(config, {MathJax: {Hub: {Config: value => settings = value}}});
      assert.equal(settings.extensions[1], root + 'static/xypic.js');
      const xypic = fs.readFileSync(publicDir + '/static/xypic.js', 'utf8');
      let completed;
      vm.runInNewContext(xypic.slice(xypic.lastIndexOf('\nMathJax.Ajax.loadComplete')), {
        document: {documentElement: {getAttribute: () => root}},
        MathJax: {Ajax: {loadComplete: value => completed = value}},
      });
      assert.equal(completed, settings.extensions[1]);
      const templates = html.match(/<script type="application\/json" id="type-templates">([\s\S]*?)<\/script>/)[1];
      const code = fs.readFileSync(publicDir + '/static/dynkin-page.js', 'utf8');
      for (const route of ['A22', 'B37/', 'missing', 'C6/30']) {
        const target = {innerHTML: 'not found'};
        const document = {
          documentElement: {dataset: {siteRoot: root}},
          getElementById: id => id === 'type-templates' ? {textContent: templates} : target,
        };
        vm.runInNewContext(code, {document, location: {pathname: root + route},
          MathJax: {Hub: {Config() {}, Queue() {}}}});
        if (route === 'A22' || route === 'B37/') {
          const type = route.replace(/\/$/, '');
          assert.match(target.innerHTML, /<svg class="dynkin"/);
          assert.ok(target.innerHTML.includes('href="' + root + type + '/1"'));
          assert.ok(document.title.includes(type));
        } else {
          assert.equal(target.innerHTML, 'not found');
        }
      }
    '''
    subprocess.run(["node", "-e", script], input=json.dumps({"publicDir": str(public.resolve()), "root": site_root}),
                   text=True, check=True)
    print("Verified fallback routing, generated links, and MathJax configuration.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--public-dir", type=Path, default=ROOT / "public")
    parser.add_argument("--base-url", default="/")
    args = parser.parse_args()
    public = args.public_dir
    site_root = urlsplit(args.base_url).path.rstrip("/") + "/"
    count, longest, large = 0, 0, 0
    models = list(records())
    expected_pages = render_pages(models)
    model_routes = {record["route"] for record in models}
    for record in models:
        page = compare(record["route"], public, expected_pages[record["route"]], site_root)
        assert json.loads(page.model) == exact_json(record), record["route"]
        for value in record["invariants"].values():
            for number in value if isinstance(value, list) else [value]:
                if isinstance(number, int):
                    longest = max(longest, len(str(number)))
                    large += abs(number) > 2**53 - 1
        count += 1
    for route, expected_html in expected_pages.items():
        if route not in model_routes:
            compare(route, public, expected_html, site_root)
            count += 1
    for asset in (ROOT / "static").rglob("*"):
        if asset.is_file():
            path = asset.relative_to(ROOT)
            assert (public / path).read_bytes() == asset.read_bytes(), (path, "asset differs")
    for alias in ("B1/1", "C1/1"):
        source = (public / alias / "index.html").read_text()
        assert 'http-equiv="refresh"' in source and "/A1/1/" in source
    certificates = "references/quantum-certificates.json"
    assert (public / certificates).read_bytes() == (Path(__file__).parents[1] / certificates).read_bytes()
    vogel = (public / "vogel/index.html").read_text()
    point_links = set(re.findall(r'xlink:href="(/[^"#]+)"', vogel))
    assert point_links
    for route in point_links:
        assert route.startswith(site_root), route
        assert (public / route[len(site_root):] / "index.html").is_file(), (route, "missing Vogel point target")
    fallback = (public / "404.html").read_text()
    check_site_links(fallback, site_root)
    assert f'data-site-root="{site_root}"' in fallback
    assert f'src="{site_root}static/dynkin-page.js"' in fallback
    verify_fallback(public, site_root)
    print(f"Verified {count} complete page bodies against Jinja content, including diagrams, controls, scripts, and citations.")
    print(f"Checked model JSON, assets, aliases, {len(point_links)} Vogel point links, and the quantum-certificate download.")
    print(f"Longest integer: {longest} digits; {large} invariant values exceed 2**53 - 1.")
    verify_dynamic_types()


if __name__ == "__main__":
    main()
