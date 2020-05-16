import json

R.<t> = PolynomialRing(QQ)

def index(D, k):
    R = RootSystem(D)
    L = R.root_lattice()

    beta = L.simple_root(k)

    length = beta.to_ambient().dot_product(beta.to_ambient())
    betabar = [alpha for alpha in L.roots()
        if alpha.coefficient(k) == 1
        and alpha.to_ambient().dot_product(alpha.to_ambient()) == length][-1]

    M = R.ambient_space()
    rho = M.rho()
    xi = M.fundamental_weight(k)
    beta = beta.to_ambient()
    betabar = betabar.to_ambient()

    return rho.dot_product(beta + betabar) / xi.dot_product(beta)


def rank(D, k):
    if D == "A1": return 2

    D = DynkinDiagram(D)
    T = [r for r in D.vertices() if r != k]

    WG = D.root_system().root_lattice().weyl_group()
    WL = D.subtype(T).root_system().root_lattice().weyl_group()

    return WG.cardinality() / WL.cardinality()


# give the "canonical" quotient, where G is the Aut^0
def canonical(D, k):
    T = D[0]
    n = int(D[1:])

    # the exceptional identifications
    if T == "B" and k == n and n >= 3: return ("D" + str(n + 1), n + 1)
    if T == "C" and k == 1: return ("A" + str(2*n - 1), 1)
    if T == "G" and k == 1: return ("B3", 1)

    return (D, k)


def dimension(D, k=0):
    def dimGB(D):
        if D == "A1": return (3, 2)

        dimG = len(D.root_system().root_poset()) * 2 + D.rank()
        dimB = len(D.root_system().root_poset()) + D.rank()
        return (dimG, dimB)

    D = DynkinDiagram(D)

    # we care about the algebraic group only
    if k == 0:
        return dimGB(D)[0]
    # we care about the quotient
    else:
        if D == DynkinDiagram("A1"): return 1

        I = [i for i in D.vertices() if i != k]

        return dimGB(D)[0] - dimGB(D)[1] - (dimGB(D.subtype(I))[0] - dimGB(D.subtype(I))[1])


def betti(D, k):
    if D == "A1": return [1, 1]

    R.<x> = PolynomialRing(ZZ)

    D = DynkinDiagram(D)
    I = [i for i in D.vertices() if i != k]
    WG = D.root_system().root_lattice().weyl_group()
    WL = D.subtype(I).root_system().root_lattice().weyl_group()

    P = prod([sum([x^i for i in range(n)]) for n in WG.degrees()])
    Q = prod([sum([x^i for i in range(n)]) for n in WL.degrees()])

    return R(P/Q).coefficients()


def embedding(D, k):
    return hilbert_polynomial(D, k)(1) - 1


def degree(D, k):
    return hilbert_polynomial(D, k).leading_coefficient() * factorial(dimension(D, k))


def hilbert_polynomial(D, k):
    return prod(hilbert_polynomial_factors(D, k))


def hilbert_polynomial_factors(D, k):
    M = RootSystem(D).ambient_space()
    l = M.fundamental_weights()[k] # G/P_k for now, can be replaced by arbitrary weight later

    return [1 + t * l.dot_product(a) / M.rho().dot_product(a) for a in M.positive_roots() if l.dot_product(a) != 0]


def pizero(D, k):
    T = D[0]
    n = int(D[1:])

    if T == "A":
        if n == 2*k - 1 and n >= 2: return "\\mathbb{Z}/2\\mathbb{Z}"
        else: return "1"
    elif T == "B":
        return "1"
    elif T == "C":
        return "1"
    elif T == "D":
        if n == 4 and k in [1, 3, 4]: return "\\mathbb{Z}/3\\mathbb{Z}"
        elif k in [n - 1, n]: return "1"
        else: return "\\mathbb{Z}/2\\mathbb{Z}"
    elif T == "E":
        if n == 6 and k in [2, 4]: return "\\mathbb{Z}/2\\mathbb{Z}"
        else: return "1"
    elif T == "F":
        return "1"
    elif T == "G":
        return "1"



grassmannians = []

def grassmannian(D, k):
    return {
            "type" : D,
            "parabolic" : int(k),
            "canonical" : [canonical(D, k)[0], int(canonical(D, k)[1])],

            "dimension" : int(dimension(D, k)),
            "index" : int(index(D, k)),
            "rank" : int(rank(D, k)),
            "betti": list(map(int, (betti(D, k)))),

            "embedding" : int(embedding(D, k)),
            "degree" : int(degree(D, k)),
            "hilbert_polynomial" : {
                "latex" : latex(hilbert_polynomial(D, k)),
                "plaintext" : str(hilbert_polynomial(D, k)),
                "factors" : [str(f.leading_coefficient()) for f in hilbert_polynomial_factors(D, k)],
                "factors_latex" : [latex(f.leading_coefficient()) for f in hilbert_polynomial_factors(D, k)],
                "evaluation" : [int(hilbert_polynomial(D, k)(i)) for i in range(20)]
            },

            "Aut^0" : canonical(D, k)[0],
            "pizero" : pizero(canonical(D, k)[0], canonical(D, k)[1]),
            "dim Aut" : int(dimension(canonical(D, k)[0])),
        }


# determining which G/P to consider
types = [("A1", 1)]
for prefix in ["A", "B", "C", "D"]:
    for n in range(2, 11):
        # let's not do D3 and D4
        if prefix == "D" and n < 4: continue

        D = prefix + str(n)

        for k in DynkinDiagram(D).vertices():
            types.append((D, k),)

exceptional = ["E6", "E7", "E8", "F4", "G2"]
for D in exceptional:
    for k in DynkinDiagram(D).vertices():
        types.append((D, k),)

# generate the Grassmannians
for (D, k) in types:
    G = grassmannian(D, k)
    grassmannians.append(G)

with open("grassmannians.json", "w") as f:
    json.dump(grassmannians, f, indent=4)

