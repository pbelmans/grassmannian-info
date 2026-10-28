"""Cited results for Grassmannians over an algebraically closed field of characteristic zero.

Keep theorem hypotheses here, rather than inferring results from plotted numerical spectra.
Bibliographical records live in references/catalog.json.
"""

import json
import math
from functools import lru_cache
from pathlib import Path


REFERENCES = json.loads((Path(__file__).parent / "references/catalog.json").read_text())
DUAL_DEGREES = json.loads((Path(__file__).parent / "references/projective-duality.json").read_text())


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
        return dict(fact("Minuscule/cominuscule case, including quadric realisations.",
                         "minuscule-qh"), status="semisimple")
    if G.is_adjoint() and not G.is_coadjoint():
        return dict(fact("Adjoint, non-coadjoint case.",
                         "adjoint-qh"), status="semisimple")
    if (T, n, k) == ("C", 4, 3):
        return dict(fact("", "coadjoint-qh", "Remark 1.2"), status="semisimple")
    if T == "E" and (n, k) in {(7, 2), (7, 4), (7, 5), (8, 6)}:
        text = "Semisimple at $q=1$; "
        if (n, k) in {(7, 4), (8, 6)}:
            text += r"$c_1\star$ has distinct, nonzero eigenvalues."
        else:
            multiplicity = 2 if k == 2 else 12
            text += (r"$c_1\star$ has the eigenvalue $0$ "
                     f"with algebraic multiplicity {multiplicity}.")
        return dict(fact(text, "holmes-muratore", "Theorems 1.2 and 6.2; Remark 6.3"),
                    status="semisimple")
    obstruction = quantum_obstruction(G.index(), tuple(G.betti()))
    if obstruction:
        return dict(fact("",
                         "graded-qh", "Theorem 1.2"), status="non-semisimple", witness=obstruction)
    old_negative = (
        (T == "B" and 2 <= k < n and k % 2 == 1)
        or (T == "C" and 2 <= k < n and k % 2 == 0)
        or (T == "D" and 2 <= k < n - 1 and k % 2 == 0)
        or (T == "E" and (n, k) in {(6, 2), (6, 4), (7, 1), (7, 3), (7, 6),
                                                    (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 7), (8, 8)})
        or (T == "F" and k in {3, 4}))
    if old_negative:
        return dict(fact("",
                         "adjoint-qh", "table on p. 326"), status="non-semisimple")
    if (T, n, k) in {("B", 6, 4), ("B", 7, 4), ("B", 7, 6),
                     ("C", 5, 3), ("C", 6, 3), ("C", 6, 5),
                     ("E", 6, 3), ("E", 6, 5), ("F", 4, 2)}:
        return dict(fact("",
                         "quantum-computation"), status="semisimple")
    return {"status": "unknown", "text": ""}


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
                           "Lefschetz centre: a spinor bundle and the structure sheaf.",
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
                           "Maximal length; fullness conjectural.", "Section 5", full=False)]
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
        return [fact("HP dual: the zero category, for the full linear embedding and Beilinson decomposition.", "hpd", "Section 8")]
    if T == "quadric":
        text = ("HP dual: the classical dual quadric, for the spinor Lefschetz structure."
                if n % 2 == 0 else
                "HP dual: the double cover of the dual projective space branched over the dual quadric, for the spinor Lefschetz structure.")
        return [fact(text, "quadric-hpd", "Theorem 1.1")]
    if T == "spinor" and n == 5:
        return [fact("HP dual: the other spinor tenfold in the dual half-spin representation, for the rectangular decomposition.", "kuznetsov-hyperplanes", "Section 6.2")]
    if T == "C" and n == k == 3:
        return [fact("HP dual: the quartic itself with an Azumaya algebra, over the complement of its singular locus in the dual projective space.", "kuznetsov-hyperplanes", "Section 6.3")]
    if T == "G" and k == 2:
        return [fact("HP dual: a double cover branched along the dual sextic, with an Azumaya algebra, over an open subset of the dual projective space.", "kuznetsov-hyperplanes", "Section 6.4")]
    if T == "A":
        N = n + 1
        k = min(k, N - k)
        if k == 2 and N == 5:
            return [fact("HP dual: the Grassmannian of two-planes in the dual five-dimensional vector space.", "kuznetsov-hyperplanes", "Section 6.1")]
        if k == 2 and N in {6, 7}:
            return [fact("HP dual: a noncommutative resolution of the Pfaffian variety of skew forms of rank at most four.", "grassmannian-hpd", "Theorem 1")]
        if k == 2 and N % 2:
            return [fact("HP dual: a noncommutative resolution of the dual Pfaffian variety.", "hori-hpd", "Theorem 4.33")]
        if k == 2:
            return [fact("Conjectural HP dual: Pirozhkov's construction from pairs of a degenerate skew form and a line in its kernel.", "pirozhkov-hpd")]
        if math.gcd(k, N) == 1:
            return [fact("HP dual: matrix factorizations on a noncommutative crepant resolution, for Fonarev's rectangular decomposition.", "doyle-hpd", "Theorem 5.13 (coprime case)")]
    return []


def page_references(G):
    keys = {ref["key"] for ref in exceptional_references(G)}
    keys.update({"big-qh", "higher-fano"})
    small = small_quantum(G)
    if "reference" in small: keys.add(small["reference"])
    if G.is_coadjoint(): keys.add("coadjoint-qh")
    for result in lefschetz_collections(G) + hpd_results(G):
        keys.add(result["reference"])
    if G.D.T == "C" and not lefschetz_collections(G): keys.add("kuznetsov-polishchuk")
    return ordered_references(keys)


def projective_duality(G):
    """Cited descriptions, with exact Chern-number computations through rank 8."""
    result = _projective_duality(G)
    computed = DUAL_DEGREES.get(G.url())
    if computed and not result.get("empty"):
        if result.get("degree") is None:
            result.update(degree_computed=True, degree_reference="tevelev",
                          degree_locator="Theorems 7.1–7.2 (arXiv version)")
        result["degree"] = computed["degree"]
    return result


def _projective_duality(G):
    """Classical duality for the minimal ample-generator embedding, not an adjoint re-embedding."""
    T, n, k = presentation(G)
    result = {"reference": "tevelev", "locator": "Theorems 9.19 and 10.22 (arXiv version)"}
    if T == "projective":
        return dict(result, label="empty", description="The dual is empty: no hyperplane contains a projective tangent space of the full projective space.", empty=True)
    if T == "quadric":
        return dict(result, label="smooth quadric", description="A smooth quadric in the dual projective space.", codimension=1, degree=2, self_dual=True)
    if T == "A" and min(k, n + 1 - k) == 2:
        N = n + 1
        return dict(result, label=f"skew forms of rank at most {N - 2 if N % 2 == 0 else N - 3}", description=f"The Pfaffian locus of skew forms of rank at most {N - 2 if N % 2 == 0 else N - 3}.",
                    codimension=1 if N % 2 == 0 else 3, degree=N // 2 if N % 2 == 0 else G.degree() if N == 5 else None,
                    self_dual=N in {4, 5})
    if T == "spinor" and n == 5:
        return dict(result, label="spinor tenfold", description="The other spinor tenfold in the dual half-spin projective space.", codimension=5, degree=G.degree(), self_dual=True)
    if T == "E" and n == 6 and k in {1, 6}:
        return dict(result, label="Cartan cubic", description="The Cartan cubic hypersurface.", codimension=1, degree=3, self_dual=False,
                    degree_reference="magic-square", degree_locator="Section 4.1")
    if ((T == "C" and n == k == 3) or (T == "A" and n == 5 and k == 3)
            or (T == "spinor" and n == 6) or (T == "E" and n == 7 and k == 7)):
        group, representation = {
            "C": (r"\mathrm{Sp}_6", r"the 14-dimensional kernel of the symplectic contraction $\bigwedge^3\mathbb C^6\to\mathbb C^6$"),
            "A": (r"\mathrm{SL}_6", r"the 20-dimensional representation $\bigwedge^3\mathbb C^6$"),
            "spinor": (r"\mathrm{Spin}_{12}", r"a 32-dimensional half-spin representation of $\mathrm{Spin}_{12}$"),
            "E": (r"E_7", "the 56-dimensional minuscule representation of $E_7$"),
        }[T]
        description = (
            rf"The dual is the hypersurface $Q=0$ in $\mathbb P(V^\vee)\cong\mathbb P^{{{G.embedding()}}}$, "
            rf"where $V$ is {representation}. "
            rf"The polynomial $Q$ is the unique, up to scale, ${group}$-invariant quartic. "
            rf"The invariant symplectic form on $V$ identifies this dual with the tangential variety of ${G.latex()}$: "
            "the union of its embedded projective tangent spaces."
        )
        return dict(result, label=rf"${group}$-invariant quartic", description=description, codimension=1, degree=4, self_dual=False,
                    degree_reference="magic-square", degree_locator="Propositions 5.4, 5.5 and 5.8 (arXiv version)")
    if T == "G" and k == 2:
        return dict(result, label=r"$G_2$-invariant sextic", description="The dual sextic hypersurface.", codimension=1, degree=6, self_dual=False,
                    degree_reference="kuznetsov-hyperplanes", degree_locator="Section 6.4")
    return dict(result, label="hypersurface", description="The projective dual is a hypersurface.", codimension=1,
                degree=None, self_dual=False)
