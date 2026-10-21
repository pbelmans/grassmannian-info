"""Check the Python -> JSON -> Hugo -> HTML boundary for every exported page."""

from html.parser import HTMLParser
import json
from pathlib import Path

from export import exact_json, records


class Page(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.model = ""
        self.invariants = {}
        self.betti = []
        self.ids, self.citations = set(), set()
        self.current = None
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        if attrs.get("href", "").startswith("#reference-"):
            self.citations.add(attrs["href"][1:])
        if attrs.get("id") == "model-data":
            self.current = "model"
        elif "data-invariant" in attrs:
            self.current = attrs["data-invariant"]
            self.invariants[self.current] = ""
        elif "data-betti-degree" in attrs:
            self.current = "betti"
            self.betti.append([int(attrs["data-betti-degree"]), ""])

    def handle_endtag(self, tag):
        if tag in {"script", "dd", "td"}:
            self.current = None

    def handle_data(self, value):
        if self.current == "model":
            self.model += value
        elif self.current == "betti":
            self.betti[-1][1] += value
        elif self.current is not None:
            self.invariants[self.current] += value


def main():
    public = Path(__file__).parent / "public"
    count, longest, large = 0, 0, 0
    for record in records():
        expected = exact_json(record)
        path = public / record["route"].lstrip("/") / "index.html"
        page = Page(path.read_text())
        assert json.loads(page.model) == expected, record["route"]
        assert not page.citations - page.ids, (record["route"], page.citations - page.ids)
        for field, value in expected["invariants"].items():
            if field == "Betti numbers":
                assert page.betti == [[2 * i, b] for i, b in enumerate(value)], record["route"]
                continue
            text = ", ".join(value) if isinstance(value, list) else value
            assert page.invariants[field] == text, (record["route"], field)
        for value in record["invariants"].values():
            for number in value if isinstance(value, list) else [value]:
                if isinstance(number, int):
                    longest = max(longest, len(str(number)))
                    large += abs(number) > 2**53 - 1
        count += 1
    refs = Page((public / "references/index.html").read_text())
    assert "reference-samokhin-vdk" in refs.ids
    print(f"Verified {count} pages: all model data, displayed invariants, and local citations.")
    print(f"Longest integer: {longest} digits; {large} invariant values exceed 2**53 - 1.")


if __name__ == "__main__":
    main()
