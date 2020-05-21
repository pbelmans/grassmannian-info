import json

R.<t> = PolynomialRing(QQ)

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


grassmannians = []

def grassmannian(D, k):
    return {
            "type" : D,
            "parabolic" : int(k),

            "embedding" : int(embedding(D, k)),
            "degree" : int(degree(D, k)),
            "hilbert_series" : [int(hilbert_polynomial(D, k)(i)) for i in range(20)],
        }


# determining which G/P to consider
types = [("A1", 1)]
for prefix in ["A", "B", "C", "D"]:
    for n in range(2, 21):
        # let's not do D2
        if prefix == "D" and n < 3: continue

        D = prefix + str(n)

        for k in DynkinDiagram(D).vertices():
            types.append((D, k),)

exceptional = ["E6", "E7", "E8", "F4", "G2"]
for D in exceptional:
    for k in DynkinDiagram(D).vertices():
        types.append((D, k),)

# generate the Grassmannians
for (D, k) in types:
    print(D, k)
    G = grassmannian(D, k)
    grassmannians.append(G)

with open("grassmannians.json", "w") as f:
    json.dump(grassmannians, f, indent=4)

