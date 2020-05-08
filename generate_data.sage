import json

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
    W = WeylGroup(D)
    I = [i for i in DynkinDiagram(D).vertices() if i != k]

    cosets = set([w.coset_representative(index_set=I) for w in W])
    lengths = [len(w.reduced_word()) for w in cosets]
    d = max(lengths)

    diagonal = [0]*(d + 1)

    for i in range(d + 1):
        diagonal[i] = len([l for l in lengths if l == i])

    return diagonal


grassmannians = []

# needs special attention, script isn't 100% robust
G = {"type" : "A1", "parabolic" : 1, "dimension" : 1 , "index" : 2, "rank" : 2, "betti" : [1, 1]}
grassmannians.append(G)


def grassmannian(D, k):
    return {"type" : D, "parabolic" : k, "dimension" : dimension(D, k), "index" : index(D, k), "rank" : rank(D, k), "betti": betti(D, k)}


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
    print(D, k)
    G = grassmannian(D, k)
    print(G)
    grassmannians.append(G)

# JSON module needs actual int's
for G in grassmannians:
    for key in ["parabolic", "dimension", "index", "rank"]:
        G[key] = int(G[key])

with open("grassmannians.json", "w") as f:
    json.dump(grassmannians, f, indent=4)

