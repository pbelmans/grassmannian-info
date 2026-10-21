"""Cited results for Grassmannians over an algebraically closed field of characteristic zero.

Keep theorem hypotheses here, rather than inferring results from plotted numerical spectra.
Bibliographical records live in references/catalog.json.
"""

import json
import math
from functools import lru_cache
from pathlib import Path


REFERENCES = json.loads((Path(__file__).parent / "references/catalog.json").read_text())


def ordered_references(keys):
    return sorted((dict(REFERENCES[key], key=key) for key in set(keys)),
                  key=lambda ref: (ref["year"], ref["label"], ref["key"]))


def fact(text, reference, locator=""):
    return {"text": text, "reference": reference, "locator": locator}


@lru_cache(maxsize=None)
def quantum_obstruction(index, betti):
    """A witness to Galkin–Leung–Li–Xiong, Theorem 1.2 (arXiv v3).

    betti[j] = b_{2j}. The grading is complex degree modulo the Fano index.
    Failure proves non-semisimplicity; passing the test proves nothing.
    """
    dimensions = [sum(betti[i::index]) for i in range(index)]
    for i in range(index):
        for d in range(index):
            j = d * i % index
            if dimensions[i] > dimensions[j]:
                return {"i": i, "d": d, "j": j, "left": dimensions[i], "right": dimensions[j]}
    return None


def small_quantum(G):
    T, n, k = G.D.T, G.D.n, G.k
    if (G.is_minuscule() or G.is_cominuscule()
            or (T == "G" and k == 1)):
        return dict(fact("The minuscule/cominuscule theorem applies, also through the quadric realisation.",
                         "minuscule-qh"), status="semisimple")
    if G.is_adjoint() and not G.is_coadjoint():
        return dict(fact("This is an adjoint, non-coadjoint variety.",
                         "adjoint-qh"), status="semisimple")
    obstruction = quantum_obstruction(G.index(), tuple(G.betti()))
    if obstruction:
        return dict(fact("The index and Betti numbers violate the necessary graded-dimension inequality.",
                         "graded-qh", "Theorem 1.2"), status="non-semisimple", witness=obstruction)
    old_negative = (
        (T == "B" and 2 <= k < n and k % 2 == 1)
        or (T == "C" and 2 <= k < n and k % 2 == 0)
        or (T == "D" and 2 <= k < n - 1 and k % 2 == 0)
        or (T == "E" and (n, k) in {(6, 2), (6, 4), (7, 1), (7, 3), (7, 6),
                                                    (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 7), (8, 8)})
        or (T == "F" and k in {3, 4}))
    if old_negative:
        return dict(fact("Non-semisimplicity is recorded in the classification table.",
                         "adjoint-qh", "table on p. 326"), status="non-semisimple")
    if (T, n, k) in {("B", 6, 4), ("B", 7, 4), ("B", 7, 6),
                     ("C", 5, 3), ("C", 6, 3), ("C", 6, 5),
                     ("E", 6, 3), ("E", 6, 5), ("F", 4, 2)}:
        return dict(fact("Quantum multiplication by the hyperplane class has square-free characteristic polynomial at q = 1; see the exact computation linked below.",
                         "quantum-computation"), status="semisimple")
    return {"status": "unknown", "text": "No conclusion is recorded here. Passing the graded-dimension test does not establish semisimplicity."}


def exceptional_references(G):
    refs = {"samokhin-vdk"}
    for D, k in G.isomorphisms():
        T, n = D[0], int(D[1:])
        if T == "A":
            refs.add("kapranov")
            if k in {1, n}: refs.add("beilinson")
        if T in {"B", "D"} and k == 1: refs.add("kapranov")
        if T == "B" and k == 2: refs.add("kuznetsov-lines")
        if T == "C":
            refs.add("guseva-novikov")
            if k == 2: refs.add("kuznetsov-lines")
            if k == n: refs.add("fonarev-lagrangian")
            if (n, k) == (3, 3): refs.add("samokhin-lg3")
            if (n, k) == (4, 3): refs.add("guseva-38")
            if k == n and n in {4, 5}: refs.add("polishchuk-samokhin")
            if k == 3: refs.add("novikov-3")
        if T == "D":
            if k == 2 and n >= 4: refs.add("kuznetsov-smirnov-d")
            if n == 5 and k in {4, 5}: refs.add("kuznetsov-hyperplanes")
            if n == 6 and k in {5, 6}: refs.add("spinor15")
        if T == "E" and n == 6 and k in {1, 6}: refs.add("faenzi-manivel")
        if T == "F":
            if k == 1: refs.add("smirnov-f4")
            if k == 4: refs.add("belmans-kuznetsov-smirnov")
        if T == "G" and k == 2: refs.add("kuznetsov-hyperplanes")
    return ordered_references(refs)


def presentation(G):
    """Recognise low-rank isomorphisms without changing the chosen polarisation."""
    forms = sorted(G.isomorphisms())
    for D, k in forms:
        n = int(D[1:])
        if D[0] == "A" and k in {1, n}:
            return "projective", n, 1
    for D, k in forms:
        n = int(D[1:])
        if D[0] in {"B", "D"} and k == 1:
            return "quadric", G.dimension(), 1
        if D == "C2" and k == 2:
            return "quadric", 3, 1
    for D, k in forms:
        n = int(D[1:])
        if D[0] == "D" and k in {n - 1, n}:
            return "spinor", n, n
    return G.D.T, G.D.n, G.k


def collection(support, reference, construction, locator="", minimal=None, full=True):
    return dict(fact(construction, reference, locator), support=support, full=full,
                rectangular=len(set(support)) == 1, minimal=minimal)


def lefschetz_collections(G):
    T, n, k = presentation(G)
    if T == "projective":
        return [collection([1] * (n + 1), "beilinson", r"The line bundles $\mathcal O,\ldots,\mathcal O(n)$.", minimal=True)]
    if T == "quadric":
        support = [2] + [1] * (n - 1) if n % 2 else [2, 2] + [1] * (n - 2)
        return [collection(support, "quadric-hpd",
                           "A chosen spinor bundle and the structure sheaf form the Lefschetz centre; this is the structure used for the HP dual below.",
                           "Lemma 2.4", minimal=True)]
    if T == "A":
        N = n + 1
        k = min(k, N - k)
        if k == 2:
            m = N // 2
            support = [m] * N if N % 2 else [m] * m + [m - 1] * m
            return [collection(support, "kuznetsov-lines", r"Symmetric powers of the dual tautological bundle form the first block.", "Theorem 4.1", minimal=True)]
        if math.gcd(k, N) == 1:
            return [collection([math.comb(N, k) // N] * N, "fonarev-minimal",
                               "Schur functors indexed by upper-triangular Young diagrams.", "Theorem 4.1 and Proposition 4.8", minimal=True)]
        # Group Kapranov's partitions in the k × (N-k) rectangle by their last row.
        return [collection([math.comb(N - t - 1, k - 1) for t in range(N - k + 1)],
                           "kapranov", r"Kapranov's Schur bundles, grouped by the last row of the Young diagram (determinant twists).",
                           minimal=False)]
    if T == "B" and k == 2:
        return [collection([n] * (2 * n - 2), "kuznetsov-lines",
                           "Symmetric powers of the dual tautological bundle and a spinor bundle.", "Theorem 7.1", minimal=True)]
    if T == "C":
        if k == 2:
            return [collection([n] * (n - 1) + [n - 1] * n, "kuznetsov-lines",
                               "Symmetric powers of the dual tautological bundle; the highest power occurs in the first n − 1 blocks.", "Theorem 5.1", minimal=True)]
        if k == 3:
            width = 2 * n - 2
            if n % 3 == 1:
                m = (n - 1) // 3
                height = (G.rank() + 4 * m) // width
                support = [height] * (2 * m) + [height - 1] * (4 * m)
            else:
                support = [G.rank() // width] * width
            return [collection(support, "novikov-3",
                               "Equivariant bundles built from Schur functors of the tautological bundle, symplectic exterior powers, and staircase complexes.",
                               "Theorem 1 and Remark 2", minimal=True)]
        if k == n and n in {4, 5}:
            support = [4, 3, 3, 3, 3] if n == 4 else [6, 6, 5, 5, 5, 5]
            return [collection(support, "polishchuk-samokhin",
                               "Tautological bundles and equivariant extensions.", minimal=True)]
    if T == "D" and k == 2:
        return [collection([n + 3] + [n + 1] * (n - 3) + [n] * (n - 1),
                           "kuznetsov-smirnov-d", "Mutations of symmetric powers, two spinor bundles, and two additional equivariant bundles.", "Remark 3.2")]
    if T == "spinor":
        if n == 5:
            return [collection([2] * 8, "kuznetsov-hyperplanes",
                               "The structure sheaf and the dual tautological bundle, with eight twists.", "Section 6.2", minimal=True)]
        if n == 6:
            return [collection([4, 4] + [3] * 8, "spinor15",
                               "The structure sheaf, the dual tautological bundle, and two equivariant extension bundles.", "Theorem 1.1", minimal=True)]
    if T == "E" and n == 6 and k in {1, 6}:
        return [collection([3] * 3 + [2] * 9, "faenzi-manivel",
                           "Three homogeneous bundles in the first three blocks; two in the remaining blocks.", "Theorem 1", minimal=True)]
    if T == "E" and n == 7 and k == 7:
        return [collection([4, 4] + [3] * 16, "spinor15",
                           "The proposed Freudenthal collection has maximal length. Its fullness is conjectural; the general existence theorem does not prove fullness of this particular collection.", "Section 5", full=False)]
    if T == "F" and k == 1:
        return [collection([3] * 8, "smirnov-f4",
                           "The structure sheaf, a homogeneous bundle, and the extension of the tangent bundle by the structure sheaf.", "Theorem 1.1", minimal=True)]
    if T == "F" and k == 4:
        return [collection([3, 3] + [2] * 9, "belmans-kuznetsov-smirnov",
                           "Restriction of the Cayley-plane collection to a general hyperplane section.", "Theorem 1.1", minimal=True)]
    if T == "G" and k == 2:
        return [collection([2] * 3, "kuznetsov-hyperplanes",
                           "The structure sheaf and a rank-two homogeneous bundle.", "Section 6.4", minimal=True)]
    return []


def hpd_results(G):
    T, n, k = presentation(G)
    if T == "projective":
        return [fact("For the full linear embedding and the standard Beilinson decomposition the HP dual category is zero.", "hpd", "Section 8")]
    if T == "quadric":
        text = ("For the spinor Lefschetz structure, the HP dual is the classical dual quadric."
                if n % 2 == 0 else
                "For the spinor Lefschetz structure, the HP dual is the double cover of the dual projective space branched over the classical dual quadric.")
        return [fact(text, "quadric-hpd", "Theorem 1.1")]
    if T == "spinor" and n == 5:
        return [fact("The HP dual of the spinor tenfold is the other spinor tenfold in the dual half-spin representation, for its rectangular decomposition.", "kuznetsov-hyperplanes", "Section 6.2")]
    if T == "C" and n == k == 3:
        return [fact("Over the complement of the singular locus of the dual quartic, Kuznetsov's HP-dual model is the quartic itself with a sheaf of Azumaya algebras. The construction is over this open subset of the dual projective space.", "kuznetsov-hyperplanes", "Section 6.3")]
    if T == "G" and k == 2:
        return [fact("Kuznetsov constructs the HP-dual model over an open subset of the dual projective space as a double cover branched along the dual sextic, with an Azumaya algebra.", "kuznetsov-hyperplanes", "Section 6.4")]
    if T == "A":
        N = n + 1
        k = min(k, N - k)
        if k == 2 and N == 5:
            return [fact("The HP dual is the Grassmannian of two-planes in the dual five-dimensional vector space.", "kuznetsov-hyperplanes", "Section 6.1")]
        if k == 2 and N in {6, 7}:
            return [fact("The HP dual is a noncommutative resolution of the Pfaffian variety of skew forms of rank at most four.", "grassmannian-hpd", "Theorem 1")]
        if k == 2 and N % 2:
            return [fact("Rennemo–Segal prove Pfaffian–Grassmannian HPD via a noncommutative resolution of the dual Pfaffian variety.", "hori-hpd", "Theorem 4.33")]
        if k == 2:
            return [fact("Pirozhkov constructs a candidate using pairs consisting of a degenerate skew form and a line in its kernel. The identification with the HP-dual category is conjectural.", "pirozhkov-hpd")]
        if math.gcd(k, N) == 1:
            return [fact("For the rectangular Fonarev decomposition, Doyle constructs the HP dual as a category of matrix factorizations on a noncommutative crepant resolution.", "doyle-hpd", "Theorem 5.13 (coprime case)")]
    return []


def page_references(G):
    keys = {ref["key"] for ref in exceptional_references(G)}
    keys.update({"big-qh", "higher-fano"})
    small = small_quantum(G)
    if "reference" in small: keys.add(small["reference"])
    if small.get("reference") == "quantum-computation": keys.add("quantum-chevalley")
    if G.is_coadjoint(): keys.add("coadjoint-qh")
    for result in lefschetz_collections(G) + hpd_results(G):
        keys.add(result["reference"])
    if G.D.T == "C" and not lefschetz_collections(G): keys.add("kuznetsov-polishchuk")
    return ordered_references(keys)
