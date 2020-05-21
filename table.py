from flask import Flask
from flask import render_template, redirect

import json
import math
import numpy as np

app = Flask(__name__)


class Dynkin:
    def __init__(self, T, n):
        assert T in ["A", "B", "C", "D", "E", "F", "G"]
        # and more consistency checks if you like

        self.T = T
        self.n = n

    def cartan_matrix(self):
        M = [[2 if i == j else 0 for i in range(self.n)] for j in range(self.n)]

        if self.T == "A":
            for i in range(self.n - 1): M[i][i + 1] = M[i + 1][i] = -1
        if self.T == "B":
            for i in range(self.n - 1): M[i][i + 1] = M[i + 1][i] = -1
            M[self.n - 2][self.n - 1] = -2
        if self.T == "C":
            for i in range(self.n - 1): M[i][i + 1] = M[i + 1][i] = -1
            M[self.n - 1][self.n - 2] = -2
        if self.T == "D":
            for i in range(self.n - 2): M[i][i + 1] = M[i + 1][i] = -1
            M[self.n - 3][self.n - 1] = M[self.n - 1][self.n - 3] = -1
        if self.T == "E":
            for i in [0, 1]: M[i][i + 2] = M[i + 2][i] = -1
            for i in range(2, self.n - 1): M[i][i + 1] = M[i + 1][i] = -1
        if self.T == "F":
            for i in range(3): M[i][i + 1] = M[i + 1][i] = -1
            M[1][2] = -2
        if self.T == "G":
            M[0][1] = -1
            M[1][0] = -3

        return M


    def coxeter_number(self):
        if self.T == "A": return self.n + 1
        if self.T == "B": return 2 * self.n
        if self.T == "C": return 2 * self.n
        if self.T == "D": return 2 * self.n - 2
        if self.T == "E" and self.n == 6: return 12
        if self.T == "E" and self.n == 7: return 18
        if self.T == "E" and self.n == 8: return 30
        if self.T == "F": return 12
        if self.T == "G": return 6


    def exponents(self):
        if self.T == "A": return list(range(1, self.n + 1))
        if self.T == "B": return list(range(1, 2*self.n, 2))
        if self.T == "C": return list(range(1, 2*self.n, 2))
        if self.T == "D": return sorted(list(range(1, 2*self.n - 2, 2)) + [self.n - 1])
        if self.T == "E" and self.n == 6: return [1, 4, 5, 7, 8, 11]
        if self.T == "E" and self.n == 7: return [1, 5, 7, 9, 11, 13, 17]
        if self.T == "E" and self.n == 8: return [1, 7, 11, 13, 17, 19, 23, 29]
        if self.T == "F": return [1, 5, 7, 11]
        if self.T == "G": return [1, 5]


    def weyl_group_cardinality(self):
        if self.T == "A": return math.factorial(self.n + 1)
        if self.T == "B": return 2**self.n * math.factorial(self.n)
        if self.T == "C": return 2**self.n * math.factorial(self.n)
        if self.T == "D": return 2**(self.n-1) * math.factorial(self.n)
        if self.T == "E" and self.n == 6: return 51840
        if self.T == "E" and self.n == 7: return 2903040
        if self.T == "E" and self.n == 8: return 696729600
        if self.T == "F": return 1152
        if self.T == "G": return 12


    def number_of_roots(self):
        if self.T == "A": return self.n * (self.n + 1)
        if self.T == "B": return 2 * self.n**2
        if self.T == "C": return 2 * self.n**2
        if self.T == "D": return 2 * self.n * (self.n - 1)
        if self.T == "E" and self.n == 6: return 72
        if self.T == "E" and self.n == 7: return 126
        if self.T == "E" and self.n == 8: return 240
        if self.T == "F": return 48
        if self.T == "G": return 12


    def number_of_positive_roots(self):
        return self.number_of_roots() // 2


    def dimension(self):
        # dimension of the associated simple Lie algebra (and algebraic group)
        return self.number_of_roots() + self.n


    def subtype(self, k):
        assert k in range(1, self.n + 1)

        if self.T == "A":
            if k == 1 or k == self.n: return Dynkin("A", self.n - 1)
            return [Dynkin("A", k - 1), Dynkin("A", self.n - k)]
        if self.T == "B":
            if k == 1: return Dynkin("B", self.n - 1)
            if k == self.n: return Dynkin("A", self.n - 1)
            return [Dynkin("A", k - 1), Dynkin("B", self.n - k)]
        if self.T == "C":
            if k == 1: return Dynkin("C", self.n - 1)
            if k == self.n: return Dynkin("A", self.n - 1)
            return [Dynkin("A", k - 1), Dynkin("C", self.n - k)]
        if self.T == "D":
            if k == 1: return Dynkin("D", self.n - 1)
            if k == self.n - 2: return [Dynkin("A", self.n - 3), Dynkin("A", 1), Dynkin("A", 1)]
            if k in [n - 1, self.n]: return Dynkin("A", self.n - 1)
            return [Dynkin("A", k - 1), Dynkin("D", self.n - k)]
        if self.T == "E":
            if k == 1: return Dynkin("D", self.n - 1)
            if k == 2: return Dynkin("A", self.n - 1)
            if k == 3: return [Dynkin("A", 1), Dynkin("A", self.n - 2)]
            if k == 4: return [Dynkin("A", 2), Dynkin("A", 1), Dynkin("A", self.n - 4)]
            if k == 5: return [Dynkin("A", 4), Dynkin("A", self.n - 5)]
            if k == 6:
                if self.n == 6: return Dynkin("D", 5)
                return [Dynkin("D", 5), Dynkin("A", self.n - 6)]
            if k == 7:
                if self.n == 7: return Dynkin("E", 6)
                return [Dynkin("E", 6), Dynkin("A", 1)]
            return Dynkin("E", 7)
        if self.T == "F":
            if k == 1: return Dynkin("C", 3)
            if k == 4: return Dynkin("B", 3)
            return [Dynkin("A", 2), Dynkin("A", 1)]
        if self.T == "G":
            return [Dynkin("A", 1)]


    def lie(self, unicode=True):
        # textual representation of associated simple Lie algebra (Unicode or TeX)
        if self.T == "A":
            return "𝖘𝖑<sub>{}</sub>".format(self.n + 1) if unicode else "\\mathfrak{sl}_{{ {} }}".format(self.n + 1)
        if self.T == "B":
            return "𝖘𝖔<sub>{}</sub>".format(2*self.n + 1) if unicode else "\\mathfrak{so}_{{ {} }}".format(2*self.n + 1)
        if self.T == "C":
            return "𝖘𝖕<sub>{}</sub>".format(2*self.n) if unicode else "\\mathfrak{sp}_{{ {} }}".format(2*self.n)
        if self.T == "D":
            return "𝖘𝖔<sub>{}</sub>".format(2*self.n) if unicode else "\\mathfrak{so}_{{ {} }}".format(2*self.n)
        if self.T == "E":
            return "𝖊<sub>{}</sub>".format(self.n) if unicode else "\\mathfrak{e}_{{ {} }}".format(self.n)
        if self.T == "F":
            return "𝖋<sub>{}</sub>".format(self.n) if unicode else "\\mathfrak{f}_{{ {} }}".format(self.n)
        if self.T == "G":
            return "𝖌<sub>{}</sub>".format(self.n) if unicode else "\\mathfrak{g}_{{ {} }}".format(self.n)



class Grassmannian:
    def __init__(self, T, n, k):
        assert k in range(1, n + 1)

        self.D = Dynkin(T, n)
        self.k = k

        self.L = self.D.subtype(k)
        # always treat the Levi as a list of simples
        if isinstance(self.L, Dynkin):
            self.L = [self.L]


    def url(self):
        return "{}{}/{}".format(self.D.T, self.D.n, self.k)


    def dimension(self):
        result = self.D.number_of_positive_roots()
        for D in self.L:
            result = result - D.number_of_positive_roots()
        return result


    def rank(self):
        result = self.D.weyl_group_cardinality()
        for D in self.L:
            result = result // D.weyl_group_cardinality()
        return result


    def canonical(self):
        # the "canonical" quotient, such that G is Aut^0 (need a better name for this)
        if self.D.T == "B" and self.k == self.D.n and self.D.n >= 3: return Grassmannian("D", self.D.n + 1, self.D.n + 1)
        if self.D.T == "C" and self.k == 1: return Grassmannian("A", 2 * self.D.n - 1, 1)
        if self.D.T == "G" and self.k == 1: return Grassmannian("B", 3, 1)

        return self


    def pizero(self):
        # the \pi_0 of the Aut
        if self.D.T == "A":
            if self.D.n == 2*self.k - 1 and self.D.n >= 2: return "\\mathbb{Z}/2\\mathbb{Z}"
        if self.D.T == "D":
            if self.D.n == 4 and self.k in [1, 3, 4]: return "\\mathbb{Z}/3\\mathbb{Z}"
            if self.k in [self.D.n - 1, self.D.n]: return "1"
            return "\\mathbb{Z}/2\\mathbb{Z}"
        if self.D.T == "E" and self.D.n == 6 and self.k in [2, 4]: return "\\mathbb{Z}/2\\mathbb{Z}"
        return "1"


    def betti(self):
        #if self.D.T == "A" and self.D.n == 1: return [1, 1]

        P = np.poly1d([1])
        for e in self.D.exponents():
            factor = np.poly1d([0])
            for i in range(e + 1): factor = np.polyadd(factor, np.poly1d([1] + [0]*i))
            P = np.polymul(P, factor)

        Q = np.poly1d([1])
        for D in self.L:
            for e in D.exponents():
                factor = np.poly1d([0])
                for i in range(e + 1): factor = np.polyadd(factor, np.poly1d([1] + [0]*i))
                Q = np.polymul(Q, factor)

        return list(map(int, np.polydiv(P, Q)[0]))


    def index(self):
        # taken from Akhiezer, Lie group actions in complex geometry, pages 124--125
        if self.D.T == "A":
            return self.D.n + 1
        if self.D.T == "B":
            if self.k == self.D.n: return 2*self.D.n
            return 2*self.D.n - self.k
        if self.D.T == "C":
            return 2*self.D.n - self.k + 1
        if self.D.T == "D":
            if self.k in [self.D.n - 1, self.D.n]: return 2*self.D.n - 2
            return 2*self.D.n - self.k - 1
        if self.D.T == "E" and self.D.n == 6: return [12, 11, 9, 7, 9, 12][self.k - 1]
        if self.D.T == "E" and self.D.n == 7: return [17, 14, 11, 8, 10, 13, 18][self.k - 1]
        if self.D.T == "E" and self.D.n == 8: return [23, 17, 13, 9, 11, 14, 19, 29][self.k - 1]
        if self.D.T == "F": return [8, 5, 7, 11][self.k - 1]
        if self.D.T == "G": return [5, 3][self.k - 1]


    def is_cominuscule(self):
        if self.D.T == "A": return True
        if self.D.T == "B" and self.k == 1: return True
        if self.D.T == "C" and self.k == self.D.n: return True
        if self.D.T == "D" and self.k in [1, self.D.n - 1, self.D.n]: return True
        if self.D.T == "E" and self.D.n == 6 and self.k in [1, 6]: return True
        if self.D.T == "E" and self.D.n == 7 and self.k == 7: return True

        return False

    def is_minuscule(self):
        if self.D.T == "A": return True
        if self.D.T == "B" and self.k == self.D.n: return True
        if self.D.T == "C" and self.k == 1: return True
        if self.D.T == "D" and self.k in [1, self.D.n - 1, self.D.n]: return True
        if self.D.T == "E" and self.D.n == 6 and self.k in [1, 6]: return True
        if self.D.T == "E" and self.D.n == 7 and self.k == 7: return True

        return False

    def is_coadjoint(self):
        if self.D.T == "B" and self.k == 1: return True
        if self.D.T == "C" and self.k == 2: return True
        if self.D.T == "D" and self.D.n >= 4 and self.k == 2: return True
        if self.D.T == "E" and self.D.n == 6 and self.k == 2: return True
        if self.D.T == "E" and self.D.n == 7 and self.k == 1: return True
        if self.D.T == "E" and self.D.n == 8 and self.k == 8: return True
        if self.D.T == "F" and self.k == 4: return True
        if self.D.T == "G" and self.k == 1: return True

        return False

    def is_adjoint(self):
        if self.D.T == "B" and self.k == 2: return True
        if self.D.T == "C" and self.k == 1: return True
        if self.D.T == "D" and self.D.n >= 4 and self.k == 2: return True
        if self.D.T == "E" and self.D.n == 6 and self.k == 2: return True
        if self.D.T == "E" and self.D.n == 7 and self.k == 1: return True
        if self.D.T == "E" and self.D.n == 8 and self.k == 8: return True
        if self.D.T == "F" and self.k == 1: return True
        if self.D.T == "G" and self.k == 2: return True

        return False


    def latex(self):
        # TODO eventually the latex function just come here?
        return latex(self.D.T, self.D.n, self.k)

    def html(self):
        # TODO eventually the html function just come here?
        return html(self.D.T, self.D.n, self.k)

    def name(self):
        # TODO eventually the name function just come here?
        return name(self.D.T, self.D.n, self.k)


# registering the Dynkin and Grassmannian class so that Jinja2 can use it directly
app.add_template_global(Dynkin, "Dynkin")
app.add_template_global(Grassmannian, "Grassmannian")


grassmannians = {letter : {} for letter in "ABCDEFG"}

G = Grassmannian("A", 18, 10)
print(G.betti())
print(G.rank(), sum(G.betti()), len(G.betti()))
# TODO if G.rank() != sum(G.betti()) there is a computational error

def latex(T, n, k):
    if T == "A":
        if k == 1: return "\\mathbb{{P}}^{{{}}}".format(n)
        elif k == n: return "\\mathbb{{P}}^{{{},\\vee}}".format(n)
        else: return "\\Gr({},{})".format(k, n + 1)

    elif T == "B":
        if k == 1: return "\\mathrm{{Q}}^{{{}}}".format(2*n - 1)
        elif n == 2 and k == 2: return latex("A", 3, 1)
        else: return "\\OGr({},{})".format(k, 2*n + 1)

    elif T == "C":
        if k == 1: return latex("A", 2*n - 1, 1)
        elif k == n:
            if n == 2: return "\\mathrm{{Q}}^{{3}}"
            else: return "\\LGr({},{})".format(k, 2*n)
        else: return "\\SGr({},{})".format(k, 2*n)

    elif T == "D":
        if k == 1: return "\\mathrm{{Q}}^{{{}}}".format(2*n - 2)
        elif n == 3 and k in [2, 3]: return latex("A", 3, 1)
        elif n == 4 and k in [3, 4]: return latex("D", 4, 1)
        elif k <= n - 2: return "\\OGr({},{})".format(k, 2*n)
        else: return "\\OGr_+({},{})".format(n, 2*n)

    elif T == "E":
        if n == 6 and k == 1: return "\\mathbb{OP}^2"
        if n == 6 and k == 6: return "\\mathbb{OP}^{2,\\vee}"

    elif T == "G":
        if k == 1: return latex("B", 3, 1)

    return "\\mathrm{{{}}}_{{{}}}/\\mathrm{{P}}_{{{}}}".format(T, n, k)


def html(T, n, k):
    if T == "A":
        if k == 1: return "&#x2119;<sup>{}</sup>".format(n)
        elif k == n: return "&#x2119;<sup>{},&#x2228;</sup>".format(n)
        else: return "Gr({},{})".format(k, n + 1)

    elif T == "B":
        if k == 1: return "Q<sup>{}</sup>".format(2*n - 1)
        elif n == 2 and k == 2: return html("A", 3, 1)
        else: return "OGr({},{})".format(k, 2*n + 1)

    elif T == "C":
        if k == 1: return html("A", 2*n - 1, 1)
        elif k == n:
            if n == 2: return "Q<sup>3</sup>"
            else: return "LGr({},{})".format(k, 2*n)
        else: return "SGr({},{})".format(k, 2*n)

    elif T == "D":
        if k == 1: return "Q<sup>{}</sup>".format(2*n - 2)
        elif n == 3 and k in [2, 3]: return html("A", 3, 1)
        elif n == 4 and k in [3, 4]: return html("D", 4, 1)
        elif k <= n - 2: return "OGr({},{})".format(k, 2*n)
        else: return "OGr<sub>+</sub>({},{})".format(n, 2*n)

    elif T == "E":
        if n == 6 and k == 1: return "&#x1D546;&#x2119;<sup>2</sup>"
        if n == 6 and k == 6: return "&#x1D546;&#x2119;<sup>2,&#x2228;</sup>"

    elif T == "G":
        if k == 1: return html("B", 3, 1)

    return "{}<sub>{}</sub>/P<sub>{}</sub>".format(T, n, k)


def plaintext(T, n, k):
    if T == "A":
        if k == 1: return "P{}".format(n)
        elif k == n: return "P{}".format(n)
        else: return "Gr({},{})".format(k, n + 1)

    elif T == "B":
        if k == 1: return "Q{}".format(2*n - 1)
        elif n == 2 and k == 2: return plaintext("A", 3, 1)
        else: return "OGr({},{})".format(k, 2*n + 1)

    elif T == "C":
        if k == 1: return plaintext("A", 2*n - 1, 1)
        elif k == n:
            if n == 2: return "Q3"
            else: return "LGr({},{})".format(k, 2*n)
        else: return "SGr({},{})".format(k, 2*n)

    elif T == "D":
        if k == 1: return "Q{}".format(2*n - 2)
        elif n == 4 and k in [3, 4]: return plaintext("D", 4, 1)
        elif k <= n - 2: return "OGr({},{})".format(k, 2*n)
        else: return "OGr+({},{})".format(n, 2*n)

    elif T == "E":
        if n == 6 and k == 1: return "OP2"
        if n == 6 and k == 6: return "OP2"

    elif T == "G":
        if k == 1: return plaintext("B", 3, 1)

    return "{}{}/P{}".format(T, n, k)


def name(T, n, k):
    if T == "A":
        if k == 1: return "projective space"
        elif k == n: return "dual projective space".format(n)
        else: return "Grassmannian"

    elif T == "B":
        if k == 1: return "quadric"
        elif n == 2 and k == 2: return name("A", 3, 1)
        else: return "orthogonal Grassmannian"

    elif T == "C":
        if k == 1: return name("A", 2*n - 1, 1)
        elif k == n:
            if n == 2: return "quadric"
            else: return "Lagrangian Grassmannian"
        else: return "symplectic Grassmannian"

    elif T == "D":
        if k == 1: return "quadric"
        elif n == 3 and k in [2, 3]: name("A", 3, 1)
        elif n == 4 and k in [3, 4]: name("D", 4, 1)
        elif k <= n - 2: return "orthogonal Grassmannian"
        else: return "orthogonal Grassmannian, spinor variety"

    elif T == "E":
        if n == 6 and k == 1: return "Cayley plane"
        if n == 6 and k == 6: return "dual Cayley plane"
        if n == 7 and k == 7: return "Freudenthal variety"

    elif T == "G":
        if k == 1: return name("B", 3, 1)
        if k == 2: return "\\mathrm{G}_2$-Grassmannian"

    return "generalised Grassmannian of type ($\\mathrm{{{}}}_{{{}}}/\mathrm{{P}}_{{{}}}$)".format(T, n, k)



# load the Sage-generated JSON file
with open("grassmannians.json") as f:
    data = json.load(f)

    for G in data:
        T = G["type"][0]
        n = int(G["type"][1:])
        k = G["parabolic"]

        # cutoff for now
        if T in "ABCD" and n > 7: continue

        G["latex"] = latex(T, n, k)
        G["html"] = html(T, n, k)
        G["plaintext"] = plaintext(T, n, k)
        G["name"] = name(T, n, k)

        # indicate type of Grassmannian
        G["cominuscule"] = G["minuscule"] = G["adjoint"] = G["coadjoint"] = False

        if T == "A": G["cominuscule"] = G["minuscule"] = True
        if T == "B" and k == 1: G["cominuscule"] = G["coadjoint"] = True
        if T == "B" and k == 2: G["adjoint"] = True
        if T == "B" and k == n: G["minuscule"] = True
        if T == "C" and k == 1: G["minuscule"] = G["adjoint"] = True
        if T == "C" and k == 2: G["coadjoint"] = True
        if T == "C" and k == n: G["cominuscule"] = True
        if T == "D" and k in [1, n - 1, n]: G["cominuscule"] = G["minuscule"] = True
        if T == "D" and k == 2: G["adjoint"] = G["coadjoint"] = True
        if T == "E" and n == 6 and k in [1, 6]: G["cominuscule"] = G["minuscule"] = True
        if T == "E" and n == 6 and k == 2: G["adjoint"] = G["coadjoint"] = True
        if T == "E" and n == 7 and k == 1: G["adjoint"] = G["coadjoint"] = True
        if T == "E" and n == 7 and k == 7: G["cominuscule"] = G["minuscule"] = True
        if T == "E" and n == 8 and k == 8: G["adjoint"] = G["coadjoint"] = True
        if T == "F" and k == 1: G["adjoint"] = True
        if T == "F" and k == 4: G["coadjoint"] = True
        if T == "G" and k == 1: G["coadjoint"] = True
        if T == "G" and k == 2: G["adjoint"] = True


        G["isomorphisms"] = [(T + str(n), k)]
        if T == "A":
            G["isomorphisms"].append(("A" + str(n), n - k + 1),)
        elif T == "B":
            if n == 2 and k == 1: G["isomorphisms"].append(("C2", 2),)
            if n == 2 and k == 2: G["isomorphisms"].extend([("C2", 1), ("A3", 1), ("A3", 3)])
            if n == 3 and k == 1: G["isomorphisms"].append(("G2", 1),)
            if n == 3 and k == 3: G["isomorphisms"].extend([("D4", 1), ("D4", 3), ("D4", 4)])
            if n == k: G["isomorphisms"].extend([("D" + str(n + 1), n), ("D" + str(n + 1), n + 1)])
        elif T == "C":
            if k == 1: G["isomorphisms"].extend([("A" + str(2*n - 1), 1), ("A" + str(2*n - 1), 2*n - 1)])
            if n == 2 and k == 1: G["isomorphisms"].append(("B2", 1),)
        elif T == "D":
            if n == 4 and k in [1, 3, 4]: G["isomorphisms"].extend([("D4", 1), ("D4", 3), ("D4", 4), ("B3", 3)])
            if k in [n - 1, n]: G["isomorphisms"].extend([("B" + str(n - 1), min(n - 1, n - (k - n + 1))), ("D" + str(n), n - (k - n + 1))])
        elif T == "E":
            if n == 6 and k in [1, 6]: G["isomorphisms"].extend([("E6", 1), ("E6", 6)])
            if n == 6 and k in [3, 5]: G["isomorphisms"].extend([("E6", 3), ("E6", 5)])
        elif T == "G":
            if k == 1: G["isomorphisms"].append(("B3", 1),)

        G["isomorphisms"] = list(set(G["isomorphisms"]))

        """
        Idea: maybe it would be better to move this to a template?

        + easier to add descriptions
        + every Lefschetz collection is a macro that can be called?

        - less convenient for MathSciNet / arXiv integration?
        - even further separation of content

        """
        G["lefschetz"] = []
        if T == "A":
            # Beilinson
            if k in [1, n]:
                G["lefschetz"].append({"support" : [1]*(n + 1)})
            # Kuznetsov (cfr. example 1.7 of http://www.mi-ras.ru/~akuznet/publications/1802.08097%5BOn%20residual%20categories%20for%20Grassmannians%5D.pdf)
            if k == 2:
                m = (n + 1) // 2
                if n % 2 == 0:
                    G["lefschetz"].append({"support" : [m]*(2*m+1)})
                else:
                    G["lefschetz"].append({"support" : [m]*m + [m-1]*m})

            # Kapranov, non-minimal

        if T == "B":
            if k == 1:
                G["lefschetz"].append({"support" : [2] + [1]*(2*n-2)})

        if T == "D":
            if k == 1 or (n == 4 and k in [3, 4]):
                G["lefschetz"].append({"support" : [2, 2] + [1]*(2*n-4)})

        if T == "E":
            # Faenzi--Manivel
            if n == 6 and k in [1, 6]:
                G["lefschetz"].append({"support" : [3]*3 + [2]*9})
        if T == "F":
            # Belmans--Kuznetsov--Smirnov
            if n == 4:
                G["lefschetz"].append({"support" : [3]*2 + [2]*8})

        # TODO check for now
        assert G["betti"] == Grassmannian(T, n, k).betti()
        assert G["index"] == Grassmannian(T, n, k).index()
        assert G["dimension"] == Grassmannian(T, n, k).dimension()

        # assigning the Grassmannian to the dictionary
        if n not in grassmannians[T]:
            grassmannians[T][n] = {}
        grassmannians[T][n][k] = G



@app.route("/")
def show_table():
    return render_template("table.html", grassmannians=grassmannians)


@app.route("/about")
def show_about():
    return render_template("about.html")


@app.route("/explained")
def show_explained():
    return render_template("explained.html")


@app.route("/freudenthal")
def show_freudenthal():
    return render_template("freudenthal.html")


@app.route("/<string:route>")
def show_type(route):
    # try for Dynkin type
    try:
        if route[0] not in ["A", "B", "C", "D", "E", "F", "G"]:
            raise ValueError

        T = route[0]
        n = int(route[1:])

        return render_template("type.html", D=Dynkin(T, n))
    except ValueError:
        pass


# TODO get rid of this one
@app.route("/<string:T><int:n>-<int:k>")
def redirect_grassmannian(T, n, k):
    return redirect("/{}{}/{}".format(T, n, k))


# TODO make these redirects?
@app.route("/<string:T><int:n>/<int:k>")
def show_grassmannian(T, n, k):
    G = grassmannians[T][n][k]

    return render_template("grassmannian.html", G=G)


@app.route("/<path:plaintext>")
def show_grassmannian_from_plaintext(plaintext):
    # look for exact plaintext match
    for T in grassmannians:
        for n in grassmannians[T]:
            for k in grassmannians[T][n]:
                G = grassmannians[T][n][k]

                if G["plaintext"] == plaintext:
                    return render_template("grassmannian.html", G=G)

    # at this point: we should look for alternative interpretations (e.g. OGr(3,7) is Q^6 etc.)
    # TODO

    # if that also fails: error
    # TODO
