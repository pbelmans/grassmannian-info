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
    D = DynkinDiagram(D)
    T = [r for r in D.vertices() if r != k]

    WG = D.root_system().root_lattice().weyl_group()
    WL = D.subtype(T).root_system().root_lattice().weyl_group()

    return WG.cardinality() / WL.cardinality()


def dimension(D, k):
    def dimGB(D):
        dimG = len(D.root_system().root_poset()) * 2 + D.rank()
        dimB = len(D.root_system().root_poset()) + D.rank()
        return (dimG, dimB)

    D = DynkinDiagram(D)
    I = [i for i in D.vertices() if i != k]

    return dimGB(D)[0] - dimGB(D)[1] - (dimGB(D.subtype(I))[0] - dimGB(D.subtype(I))[1])


def betti(D, k):
    R.<x> = PolynomialRing(ZZ)

    D = DynkinDiagram(D)
    I = [i for i in D.vertices() if i != k]
    WG = D.root_system().root_lattice().weyl_group()
    WL = D.subtype(I).root_system().root_lattice().weyl_group()

    P = prod([sum([x^i for i in range(n)]) for n in WG.degrees()])
    Q = prod([sum([x^i for i in range(n)]) for n in WL.degrees()])

    return list(map(int, R(P/Q).coefficients()))


def embedding(D, k):
    return int(hilbert_polynomial(D, k)(1))


def degree(D, k):
    return int(hilbert_polynomial(D, k).leading_coefficient() * factorial(dimension(D, k)))


def hilbert_polynomial(D, k):
    return prod(hilbert_polynomial_factors(D, k))


def hilbert_polynomial_factors(D, k):
    M = RootSystem(D).ambient_space()
    l = M.fundamental_weights()[k] # G/P_k for now, can be replaced by arbitrary weight later

    return [1 + t * l.dot_product(a) / M.rho().dot_product(a) for a in M.positive_roots() if l.dot_product(a) != 0]


grassmannians = []

# needs special attention, script isn't 100% robust
G = {
     "type" : "A1",
     "parabolic" : 1,

     "dimension" : 1,
     "index" : 2,
     "rank" : 2,
     "betti" : [int(1), int(1)],

     "embedding" : 1,
     "degree" : 1,
  }
# TODO just make the functions more robust...
grassmannians.append(G)


def grassmannian(D, k):
    print(D, k)
    return {
            "type" : D,
            "parabolic" : k,

            "dimension" : dimension(D, k),
            "index" : index(D, k),
            "rank" : rank(D, k),
            "betti": betti(D, k),

            "embedding" : embedding(D, k),
            "degree" : degree(D, k),
            "hilbert_polynomial" : {
                "latex" : latex(hilbert_polynomial(D, k)),
                "plaintext" : str(hilbert_polynomial(D, k)),
                "factors" : [str(f.leading_coefficient()) for f in hilbert_polynomial_factors(D, k)],
                "factors_latex" : [latex(f.leading_coefficient()) for f in hilbert_polynomial_factors(D, k)],
                "evaluation" : [int(hilbert_polynomial(D, k)(i)) for i in range(20)]
            }
        }


types = []
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

# JSON module needs actual int's
for G in grassmannians:
    for key in ["parabolic", "dimension", "index", "rank", "embedding", "degree"]:
        if not key in G: continue
        G[key] = int(G[key])

# debugging
#for G in grassmannians:
#    for key in G:
#        print(G, key, type(G[key]))

with open("grassmannians.json", "w") as f:
    json.dump(grassmannians, f, indent=4)

