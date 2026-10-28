"""Values revealed by the Freudenthal information toggle."""

from html.parser import HTMLParser
import unittest

from tests.support import page


class Varieties(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.invariants = []
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "td" and "data-euler" in attrs:
            self.invariants.append(tuple(int(attrs[key]) for key in
                                         ("data-dimension", "data-index", "data-euler")))


class FreudenthalTests(unittest.TestCase):
    def test_information_in_both_squares(self):
        values = Varieties(page("/freudenthal")).invariants
        self.assertEqual(len(values), 16 + 5 * 9)
        self.assertTrue(all(value > 0 for row in values for value in row))
        self.assertEqual(values[:2], [(1, 2, 2), (3, 2, 6)])
        self.assertEqual(values[4:6], [(2, 3, 3), (4, 3, 9)])
        # n=3 reproduces the upper left 3x3 square, including its custom cells.
        self.assertEqual(values[16:25], [values[i] for i in (0, 1, 2, 4, 5, 6, 8, 9, 10)])
        for n in range(3, 8):
            offset = 16 + 9 * (n - 3)
            self.assertEqual(values[offset + 3], (n - 1, n, n))
            self.assertEqual(values[offset + 4], (2 * (n - 1), n, n * n))


if __name__ == "__main__":
    unittest.main()
