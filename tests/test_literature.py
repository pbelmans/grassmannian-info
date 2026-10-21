"""Mathematical and rendering regressions for the cited classification."""
import json
from html.parser import HTMLParser
from pathlib import Path
import unittest

from application import app, data, Grassmannian
import literature


class Links(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.ids, self.citations, self.entries = set(), set(), []
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if attrs.get('href', '').startswith('#reference-'):
            self.citations.add(attrs['href'][1:])
        if 'data-exceptional' in attrs:
            self.entries.append(attrs)


class DualityTable(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.rows = {}
        self.row = None
        self.cell = None
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'tr' and 'data-duality-row' in attrs:
            self.row = attrs['id']
            self.rows[self.row] = []
        if self.row is not None and tag in {'td', 'th'}:
            self.cell = []

    def handle_data(self, value):
        if self.cell is not None:
            self.cell.append(value)

    def handle_endtag(self, tag):
        if tag in {'td', 'th'} and self.cell is not None:
            self.rows[self.row].append(''.join(self.cell).strip())
            self.cell = None
        if tag == 'tr':
            self.row = None


def varieties():
    for row in data:
        yield Grassmannian(row['type'][0], int(row['type'][1:]), row['parabolic'])


class LiteratureTests(unittest.TestCase):
    def test_full_shapes_have_correct_rank(self):
        for G in varieties():
            for c in G.lefschetz():
                with self.subTest(G=(G.D.D, G.k)):
                    self.assertEqual(sum(c['support']), G.rank())
                    self.assertEqual(sorted(c['support'], reverse=True), c['support'])
                    self.assertTrue(all(x > 0 for x in c['support']))
                    self.assertLessEqual(len(c['support']), G.index())

    def test_quantum_obstructions_do_not_contradict_positive_theorems(self):
        unknown_exceptional = set()
        for G in varieties():
            if G.small_qh_semisimple():
                self.assertIsNone(literature.quantum_obstruction(G.index(), tuple(G.betti())), (G.D.D, G.k))
            if G.D.T in 'EFG' and G.small_quantum_result()['status'] == 'unknown':
                unknown_exceptional.add((G.D.D, G.k))
        self.assertEqual(unknown_exceptional, {('E7', 2), ('E7', 4), ('E7', 5), ('E8', 6)})
        for T, n, k in [('B', 9, 8), ('C', 8, 7), ('D', 9, 7)]:
            self.assertEqual(Grassmannian(T, n, k).small_quantum_result()['reference'], 'graded-qh')
        # A passing obstruction is not a semisimplicity certificate.
        self.assertEqual(Grassmannian('E', 7, 4).small_quantum_result()['status'], 'unknown')

    def test_certificates_match_varieties(self):
        certificates = json.loads(Path('references/quantum-certificates.json').read_text())
        self.assertEqual(len(certificates), 9)
        for key, cert in certificates.items():
            D, k = key.split('/')
            G = Grassmannian(D[0], int(D[1:]), int(k))
            self.assertEqual(cert['dimension'], G.rank())
            self.assertEqual(len(cert['coefficients']), G.rank() + 1)
            self.assertEqual(cert['coefficients'][-1], '1')
            self.assertEqual(G.small_quantum_result()['reference'], 'quantum-computation')

    def test_new_collection_shapes(self):
        self.assertEqual(Grassmannian('C', 6, 3).lefschetz()[0]['support'], [16] * 10)
        self.assertEqual(Grassmannian('C', 4, 3).lefschetz()[0]['support'], [6, 6, 5, 5, 5, 5])
        self.assertEqual(Grassmannian('F', 4, 4).lefschetz()[0]['support'], [3, 3] + [2] * 9)
        self.assertFalse(Grassmannian('E', 7, 7).lefschetz()[0]['full'])
        self.assertIn('quartic itself', Grassmannian('C', 3, 3).hpd_results()[0]['text'])

    def test_front_page_is_uniform_without_result_citations(self):
        response = app.test_client().get('/')
        self.assertEqual(response.status_code, 200)
        links = Links(response.get_data(as_text=True))
        self.assertTrue(links.entries)
        self.assertTrue(all(e['data-exceptional'] == 'known' and e['data-big'] == 'semisimple' for e in links.entries))
        self.assertFalse(links.citations)

    def test_all_homogeneous_pages_and_citation_targets(self):
        client = app.test_client()
        for G in varieties():
            path = f'/{G.D.D}/{G.k}'
            with self.subTest(path=path):
                response = client.get(path)
                self.assertEqual(response.status_code, 200)
                html = response.get_data(as_text=True)
                links = Links(html)
                self.assertFalse(links.citations - links.ids)
                self.assertIn('References', html)
                self.assertNotIn('reference-tevelev', links.ids)
                self.assertNotIn('reference-magic-square', links.ids)

    def test_reference_and_legacy_routes(self):
        client = app.test_client()
        for path in ['/references', '/references/quantum-certificates.json', '/static/style.css', '/static/fanography.png', '/A5', '/horospherical/X3(3,2)']:
            with client.get(path) as response:
                self.assertEqual(response.status_code, 200, path)
        links = Links(client.get('/references').get_data(as_text=True))
        self.assertFalse(links.citations - links.ids)

    def test_classical_duality(self):
        self.assertTrue(Grassmannian('C', 6, 1).projective_duality()['empty'])
        self.assertEqual(Grassmannian('A', 6, 2).projective_duality()['codimension'], 3)
        self.assertEqual(Grassmannian('A', 5, 2).projective_duality()['degree'], 3)
        self.assertEqual(Grassmannian('B', 4, 4).projective_duality()['codimension'], 5)
        self.assertTrue(Grassmannian('D', 5, 5).projective_duality()['self_dual'])
        self.assertEqual(Grassmannian('C', 6, 3).projective_duality()['codimension'], 1)
        self.assertEqual(Grassmannian('E', 7, 7).projective_duality()['degree'], 4)

    def test_projective_duality_table_covers_rank_eight(self):
        response = app.test_client().get('/projective-duality')
        self.assertEqual(response.status_code, 200)
        html = response.get_data(as_text=True)
        rows = DualityTable(html).rows
        expected = sorted((G.D.T, G.D.n, G.k) for G in varieties() if G.D.n <= 8)
        self.assertEqual(list(rows), [f'{T}{n}-{k}' for T, n, k in expected])
        # Ambient dimension, codimension, degree and self-duality in the
        # empty, hypersurface and higher-codimension cases.
        for key, values in {
            'A3-1': ['3', '—', '—', 'no'],
            'A4-2': ['9', '3', '5', 'yes'],
            'A5-2': ['14', '1', '3', 'no'],
            'A6-2': ['20', '3', 'not recorded', 'no'],
            'C3-3': ['13', '1', '4', 'no'],
            'D5-5': ['15', '5', '12', 'yes'],
            'E7-7': ['55', '1', '4', 'no'],
        }.items():
            self.assertEqual(rows[key][2:6], values, key)
        links = Links(html)
        self.assertFalse(links.citations - links.ids)
        self.assertTrue({f'type-{T}' for T in 'ABCDEFG'} <= links.ids)


if __name__ == '__main__':
    unittest.main()
