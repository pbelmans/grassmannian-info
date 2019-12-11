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

def dimGB(D):
  dimG = len(D.root_system().root_poset()) * 2 + D.rank()
  dimB = len(D.root_system().root_poset()) + D.rank()
  return (dimG, dimB)

def dimension(D, T):
  D = DynkinDiagram(D)

  return dimGB(D)[0] - dimGB(D)[1] - (dimGB(D.subtype(T))[0] - dimGB(D.subtype(T))[1])


grassmannians = []

G = {"D" : "A1", "k" : 1, "dimension" : 1 , "index" : 2, "rank" : 2}
grassmannians.append(G)

for prefix in ["A", "B", "C", "D"]:
  for n in range(2, 11):
    if prefix == "D" and n < 4:
      continue

    D = prefix + str(n)

    for k in range(1, n + 1):
      G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, n + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
      grassmannians.append(G)

D = "E6"
for k in range(1, 7):
  G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, 6 + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
  grassmannians.append(G)

D = "E7"
for k in range(1, 8):
  G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, 7 + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
  grassmannians.append(G)

D = "E8"
for k in range(1, 9):
  G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, 8 + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
  grassmannians.append(G)

D = "F4"
for k in range(1, 5):
  G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, 4 + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
  grassmannians.append(G)

D = "G2"
for k in range(1, 3):
  G = {"D" : D, "k" : k, "dimension" : dimension(D, [i for i in range(1, 2 + 1) if i != k]), "index" : index(D, k), "rank" : rank(D, k)}
  grassmannians.append(G)

# JSON module needs actual int's
for G in grassmannians:
  for key in ["k", "dimension", "index", "rank"]:
    G[key] = int(G[key])

with open("grassmannians.json", "w") as f:
  json.dump(grassmannians, f, indent=4)

