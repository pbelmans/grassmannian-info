from flask import Flask
from flask import render_template, redirect

import os
import json
import math
import numpy as np
from datetime import datetime

app = Flask(__name__)

app.static_folder = "static"

import eigenvalues


"""
TODO
- implement UTF8 subscripts systematically
"""


class Dynkin:
    def __init__(self, T, n):
        assert T in ["A", "B", "C", "D", "E", "F", "G"]
        # and more consistency checks if you like

        self.T = T
        self.n = n
        self.D = T + str(n)

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
            if k in [self.n - 1, self.n]: return Dynkin("A", self.n - 1)
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


    def Aut(self):
        return self.canonical().D


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


    def small_qh_not_semisimple(self):
        # these are the 'no's in the table of page 33 of [MR2821244], Chaput--Perrin, On the quantum cohomology of adjoint varieties
        if self.D.T == "B" and self.k in range(2, self.D.n) and self.k % 2 == 1: return True
        if self.D.T == "C" and self.k in range(2, self.D.n) and self.k % 2 == 0: return True
        if self.D.T == "D" and self.k in range(2, self.D.n - 1) and self.k % 2 == 0: return True
        if self.D.T == "E" and self.D.n == 6 and self.k in [2, 4]: return True
        if self.D.T == "E" and self.D.n == 7 and self.k in [1, 3, 6]: return True
        if self.D.T == "E" and self.D.n == 8 and self.k in [1, 2, 3, 4, 5, 7, 8]: return True
        if self.D.T == "F" and self.k in [3, 4]: return True

        # this means that either it _is_ semisimple, or it's not known
        return False


    def small_qh_semisimple(self):
        # these are the 'yes's in the table of page 33 of [MR2821244], Chaput--Perrin, On the quantum cohomology of adjoint varieties
        if self.D.T == "A": return True
        if self.D.T == "B" and self.k in [1, 2, self.D.n]: return True # table contains a typo (?): OGr(2,2n+1) is adjoint non-coadjoint
        if self.D.T == "C" and self.k in [1, self.D.n]: return True
        if self.D.T == "D" and self.k in [1, self.D.n - 1, self.D.n]: return True
        if self.D.T == "E" and self.D.n == 6 and self.k in [1, 6]: return True
        if self.D.T == "E" and self.D.n == 7 and self.k == 7: return True
        if self.D.T == "F" and self.k == 1: return True
        if self.D.T == "G": return True

        # if all eigenvalues are distinct then we also have semisimplicity
        if self.eigenvalues():
            return len(set(self.eigenvalues())) == self.rank()

        # this means that either it is _not_ semisimple, or it's not known
        return False


    def big_qh_semisimple(self):
        # small semisimplicity implies big semisimplicity
        if self.small_qh_semisimple(): return True

        # works of Perrin, Mellit--Perrin--Smirnov, Galkin--Mellit--Smirnov
        if self.D.T == "C" and self.k == 2: return True
        # work of Perrin, see Theorem 4 in [1405.5914], Perrin: Semisimple quantum cohomology of some Fano varieties
        if self.D.T == "F" and self.k == 4: return True
        # work of Perrin--Smirnov (in progress)
        if self.is_coadjoint(): return True
        # by page 3 of 2307.03696
        if self.D.T == "C" and self.k >= 3: return True


    def latex(self):
        # TODO eventually the latex function just come here?
        return latex(self.D.T, self.D.n, self.k)

    def html(self):
        # TODO eventually the html function just come here?
        return html(self.D.T, self.D.n, self.k)

    def name(self):
        # TODO eventually the name function just come here?
        return name(self.D.T, self.D.n, self.k)

    def plaintext(self):
        # TODO eventually the plaintext function just come here?
        return plaintext(self.D.T, self.D.n, self.k)


    def isomorphisms(self):
        result = [(self.D.T + str(self.D.n), self.k)]

        if self.D.n == 1: return result

        if self.D.T == "A":
            result.append(("A" + str(self.D.n), self.D.n - self.k + 1),)
            if self.D.n == 3 and self.k in [1, 3]: result.extend([("B2", 2), ("C2", 1), ("D3", 2), ("D3", 3)])
            if self.D.n == 3 and self.k == 2: result.append(("D3", 1),)
            if self.D.n % 2 == 1 and self.k in [1, self.D.n]: result.append(("C" + str(int((self.D.n + 1) / 2)), 1),)
        if self.D.T == "B":
            if self.D.n == 2 and self.k == 1: result.append(("C2", 2),)
            if self.D.n == 2 and self.k == 2: result.extend([("D3", 2), ("D3", 3), ("C2", 1), ("A3", 1), ("A3", 3)])
            if self.D.n == 3 and self.k == 1: result.append(("G2", 1),)
            if self.D.n == 3 and self.k == 3: result.extend([("D4", 1), ("D4", 3), ("D4", 4)])
            if self.D.n == self.k: result.extend([("D" + str(self.D.n + 1), self.D.n), ("D" + str(self.D.n + 1), self.D.n + 1)])
        if self.D.T == "C":
            if self.k == 1: result.extend([("A" + str(2*self.D.n - 1), 1), ("A" + str(2*self.D.n - 1), 2*self.D.n - 1)])
            if self.D.n == 2 and self.k == 1: result.extend([("B2", 2), ("D3", 1), ("D3", 3)])
        if self.D.T == "D":
            if self.D.n == 3 and self.k in [2, 3]: result.extend([("A3", 1), ("A3", 3), ("C2", 1)])
            if self.D.n == 3 and self.k == 1: result.append(("A3", 2),)
            if self.D.n == 4 and self.k in [1, 3, 4]: result.extend([("D4", 1), ("D4", 3), ("D4", 4), ("B3", 3)])
            if self.k in [self.D.n - 1, self.D.n]: result.extend([("B" + str(self.D.n - 1), min(self.D.n - 1, self.D.n - (self.k - self.D.n + 1))), ("D" + str(self.D.n), self.D.n - (self.k - self.D.n + 1))])
        if self.D.T == "E":
            if self.D.n == 6 and self.k in [1, 6]: result.extend([("E6", 1), ("E6", 6)])
            if self.D.n == 6 and self.k in [3, 5]: result.extend([("E6", 3), ("E6", 5)])
        if self.D.T == "G":
            if self.k == 1: result.append(("B3", 1),)

        return list(set(result))

    def embedding(self):
        global grassmannians

        try:
            return grassmannians[self.D.T][self.D.n][self.k]["embedding"]
        except:
            return False

    def degree(self):
        global grassmannians

        try:
            return grassmannians[self.D.T][self.D.n][self.k]["degree"]
        except:
            return False

    def hilbert_series(self):
        global grassmannians

        try:
            return grassmannians[self.D.T][self.D.n][self.k]["hilbert_series"]
        except:
            return False

    def lefschetz(self):
        global grassmannians

        return grassmannians[self.D.T][self.D.n][self.k]["lefschetz"]

    def exceptional_sequences(self):
        return sorted(list(set([sequence for (T, k) in self.isomorphisms() for sequence in exceptional(T, k)])))#, key=lambda sequence: sequence[1])

    def eigenvalues(self):
        try:
            if self.D.T == "A": return zip(*eigenvalues.A[self.D.n, self.k])
            if self.D.T == "B": return zip(*eigenvalues.B[self.D.n, self.k])
            if self.D.T == "C": return zip(*eigenvalues.C[self.D.n, self.k])
            if self.D.T == "D": return zip(*eigenvalues.D[self.D.n, self.k])
            if self.D.T == "E": return zip(*eigenvalues.E[self.D.n, self.k])
            if self.D.T == "F": return zip(*eigenvalues.F[self.D.n, self.k])
            if self.D.T == "G": return zip(*eigenvalues.G[self.D.n, self.k])
        except KeyError: return None

        return None

    # the k-Fanoness
    def fano(self):
        # Theorems 1.3 and 1.4 of https://arxiv.org/abs/2110.02339

        # 3-Fano: Theorem 1.4 of https://arxiv.org/abs/2110.02339
        if self.D.T == "A" and self.k in [1, self.D.n] and self.D.n >= 3: return 3
        if self.D.T == "B" and self.k == 1 and self.D.n >= 4: return 3
        if self.D.T == "D" and self.k == 1 and self.D.n >= 5: return 3
        # exceptional ways of getting P^3
        if self.D.T == "B" and self.D.n == 2 and self.k == 2: return 3
        if self.D.T == "D" and self.D.n == 3 and self.k in [2, 3]: return 3
        # exceptional ways of getting P^{2n+1}
        if self.D.T == "C" and self.D.n >= 2 and self.k == 1: return 3

        # 2-Fano: Theorem 1.3 of https://arxiv.org/abs/2110.02339
        if self.D.T == "A" and self.D.n >= 2 and self.k in [1, self.D.n, math.floor((self.D.n + 1) / 2), math.ceil((self.D.n + 1) / 2)]: return 2

        if self.D.T == "B" and self.D.n >= 2 and self.k == 1: return 2
        if self.D.T == "B" and 2 * self.D.n - 1 == 3*self.k: return 2
        if self.D.T == "B" and self.D.n == self.k: return 2

        if self.D.T == "C" and self.k == 1: return 2
        if self.D.T == "C" and 2 * self.D.n + 2 == 3*self.k: return 2
        if self.D.T == "C" and self.k == self.D.n: return 2

        if self.D.T == "D" and self.k == 1: return 2
        if self.D.T == "D" and 2 * self.D.n - 2 == 3*self.k: return 2
        if self.D.T == "D" and self.k in [self.D.n - 1, self.D.n]: return 2

        if self.D.T == "E" and self.k in [1, 2, self.D.n]: return 2
        if self.D.T == "F" and self.k == 4: return 2
        if self.D.T == "G": return 2

        # just Fano, not higher Fano
        return 1


class Horospherical:
    # T and n determine the Lie type of the group G, y and z are the indices of the weights for \omega_Y and \omega_Z
    def __init__(self, T, n, y, z):
        assert T in ["B", "C", "F", "G"]
        if T == "B": assert (n == 3 and y == 1 and z == 3) or (n >= 3 and y == n - 1 and z == n)
        if T == "C": assert n >= 2 and y in range(2, n + 1) and z == y - 1
        if T == "F": assert y == 2 and z == 3
        if T == "G": assert y == 1 and z == 2

        self.D = Dynkin(T, n)
        self.y, self.z = y, z

        self.Y =  Grassmannian(T, n, y)
        self.Z =  Grassmannian(T, n, z)


    def dimension(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1: return self.D.n * (self.D.n + 3) // 2
            else: return 9
        if self.D.T == "C": return self.y * (2*self.D.n + 1 - self.y) - self.y * (self.y - 1) // 2
        if self.D.T == "F": return 23
        if self.D.T == "G": return 7


    def latex(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1: return "X^1({})".format(self.D.n)
            else: return "X^2"
        if self.D.T == "C": return "X^3({}, {})".format(self.D.n, self.y)
        if self.D.T == "F": return "X^4"
        if self.D.T == "G": return "X^5"


    def html(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1: return "<i>X</i><sup> 1</sup>({})".format(self.D.n)
            else: return "<i>X</i><sup> 2</sup>"
        if self.D.T == "C": return "<i>X</i><sup> 3</sup>({}, {})".format(self.D.n, self.y)
        if self.D.T == "F": return "<i>X</i><sup> 4</sup>"
        if self.D.T == "G": return "<i>X</i><sup> 5</sup>"


    def plaintext(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1: return "X1({})".format(self.D.n)
            else: return "X2"
        if self.D.T == "C": return "X3({},{})".format(self.D.n, self.y)
        if self.D.T == "F": return "X4"
        if self.D.T == "G": return "X5"


    def name(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1: return "horo-orthogonal Grassmannian"
            else: return "B&#x2083;-horospherical variety"
        if self.D.T == "C": return "odd symplectic Grassmannian"
        if self.D.T == "F": return "F&#x2084;-horospherical variety"
        if self.D.T == "G": return "G&#x2082;-horospherical variety"


    def codimY(self):
        return self.dimension() - self.Y.dimension()


    def codimZ(self):
        return self.dimension() - self.Z.dimension()


    def index(self):
        return self.codimY() + self.codimZ()


    def betti(self):
        Y = np.poly1d(self.Y.betti())
        Z = np.poly1d(self.Z.betti())

        PY = np.poly1d([1] * (self.codimY() + 1))
        PZ = np.poly1d([1] * (self.codimZ()))

        return np.asarray(Y * PY + Z - Z * PZ)


    def rank(self):
        return sum(self.betti())


    def exceptional_sequences(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1:
                pass
            else:
                return [("Gonzales&ndash;Pech&ndash;Perrin&ndash;Samokhin", 2018, "1803.05063")]
        if self.D.T == "C":
            if self.y == 2: return [("Pech", 2013, "MR2998953")]
            if self.D.n == 3 and self.y == 3: return [("Fonarev", 2020, "1804.06946")]
            if self.D.n == 4 and self.y == 3: return [("Cattani", 2023, "2305.06867")]
        if self.D.T == "F": pass
        if self.D.T == "G": return [("Gonzales&ndash;Pech&ndash;Perrin&ndash;Samokhin", 2018, "1803.05063")]

        return []


    def lefschetz(self):
        return []


    def eigenvalues(self):
        try:
            if self.D.T == "B":
                if self.y == self.D.n - 1:
                    return zip(*eigenvalues.X1[self.D.n])
                else:
                    return zip(*eigenvalues.X2)
            if self.D.T == "C":
                return zip(*eigenvalues.X3[self.D.n, self.y])
            if self.D.T == "F": return zip(*eigenvalues.X4)
            if self.D.T == "G": return zip(*eigenvalues.X5)
        except (AttributeError, KeyError): return None

        return None


    def small_qh_not_semisimple(self):
        return False


    def small_qh_semisimple(self):
        if self.D.T == "B":
            if self.y == self.D.n - 1:
                if self.D.n == 3:
                    return True #[("Gonzales&ndash;Pech&ndash;Perrin&ndash;Samokhin", 2018, "1803.05063")]
            else:
                return True # [("Gonzales&ndash;Pech&ndash;Perrin&ndash;Samokhin", 2018, "1803.05063")]
        if self.D.T == "C":
            if self.y == 2: return True # [("Pech", 2013, "MR2998953")]
            if self.D.n == 3 and self.y == 3: return True
        if self.D.T == "F": pass
        if self.D.T == "G": return True # [("Gonzales&ndash;Pech&ndash;Perrin&ndash;Samokhin", 2018, "1803.05063")]

        # if all eigenvalues are distinct then we also have semisimplicity
        if self.eigenvalues():
            return len(set(self.eigenvalues())) == self.rank()

        return False


    def big_qh_semisimple(self):
        # small semisimplicity implies big semisimplicity
        if self.small_qh_semisimple(): return True

    def fano(self):
        # TODO maybe not true
        return 1



# registering the Dynkin, Grassmannian and Horospherical class so that Jinja2 can use it directly
app.add_template_global(Dynkin, "Dynkin")
app.add_template_global(Grassmannian, "Grassmannian")
app.add_template_global(Horospherical, "Horospherical")


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
        elif k == n: return "P{}&vee;".format(n)
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
        elif n == 3 and k in [2, 3]: return name("A", 3, 1)
        elif n == 4 and k in [3, 4]: return name("D", 4, 1)
        elif k <= n - 2: return "orthogonal Grassmannian"
        else: return "orthogonal Grassmannian, spinor variety"

    elif T == "E":
        if n == 6 and k == 1: return "Cayley plane"
        if n == 6 and k == 6: return "dual Cayley plane"
        if n == 7 and k == 7: return "Freudenthal variety"

    elif T == "G":
        if k == 1: return name("B", 3, 1)
        if k == 2: return "$\\mathrm{G}_2$-Grassmannian"

    return "generalised Grassmannian of type {}<sub>{}</sub>/P<sub>{}</sub>".format(T, n, k)


def exceptional(D, k):
    T = D[0]
    n = int(D[1:])

    results = []
    if T == "A":
        if k in [1, n]: results.append(("Beilinson", 1978, "MR0509388"),)

        results.append(("Kapranov", 1988, "MR0939472"),)

    if T == "B":
        if k == 1: results.append(("Kapranov", 1988, "MR0939472"),)
        if k == 2: results.append(("Kuznetsov", 2008, "MR2434094"),)
        if k == 4 and n == 4: results.append(("Kuznetsov", 2006, "MR2238172"),)

    if T == "C":
        if k == 2: results.append(("Kuznetsov", 2008, "MR2434094"),)
        if k == 3 and n == 3: results.append(("Samokhin", 2001, "MR1859740"),)
        if k == 3 and n == 4: results.append(("Guseva", 2018, "1810.07777"),)
        if k == 3 and n == 5: results.append(("Novikov", 2020, ""),)
        if k == 4 and n == 4: results.append(("Polishchuk&ndash;Samokhin", 2011, "MR2822466"),)
        if k == n: results.append(("Fonarev", 2019, "1911.08968"),)

    if T == "D":
        if k == 1: results.append(("Kapranov", 1988, "MR0939472"),)
        if k == 2: results.append(("Kuznetsov&ndash;Smirnov", 2020, "2001.04148"),)
        if k == 4 and n == 5: results.append(("Kuznetsov", 2008, "MR2238172"),)
        if k == 5 and n == 6: results.append(("Benedetti&ndash;Faenzi&ndash;Smirnov", 2023, "2310.01090"))

    if T == "E":
        if k == 1 and n == 6: results.append(("Faenzi&ndash;Manivel", 2013, "MR3293722"),)

    if T == "F":
        if k == 1: results.append(("Smirnov", 2021, "2107.07814"),)
        if k == 4: results.append(("Belmans&ndash;Kuznetsov&ndash;Smirnov", 2020, "2005.01989"),)

    if T == "G":
        if k == 2: results.append(("Kuznetsov", 2008, "MR2238172"),)

    return results


grassmannians = {letter : {} for letter in "ABCDEFG"}

# load the Sage-generated JSON file
with open(os.path.join(os.path.realpath(os.path.dirname(__file__)), "grassmannians.json"), "r", encoding="utf8") as f:
    data = json.load(f)

    for G in data:
        T = G["type"][0]
        n = int(G["type"][1:])
        k = G["parabolic"]


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
            if k == 4:
                G["lefschetz"].append({"support" : [3]*2 + [2]*8})


        # assigning the Grassmannian to the dictionary
        if n not in grassmannians[T]:
            grassmannians[T][n] = {}
        grassmannians[T][n][k] = G


"""
Idea: maybe it would be better to move this to a template?

+ easier to add descriptions
+ every Lefschetz collection is a macro that can be called?

- less convenient for MathSciNet / arXiv integration?
- even further separation of content
"""


# slightly different indexing scheme: just use the plaintext X1(n), X2, X3(n,m), X4, X5
horospherical_varieties = \
        [Horospherical("B", n, n - 1, n) for n in range(3, 20)] \
        + [Horospherical("B", 3, 1, 3)] \
        + [Horospherical("C", n, m, m - 1) for n in range(2, 20) for m in range(2, n + 1)] \
        + [Horospherical("F", 4, 2, 3)] \
        + [Horospherical("G", 2, 1, 2)]
horospherical_varieties = dict({X.plaintext() : X for X in horospherical_varieties})


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

        if T in ["B", "C"] and n == 1:
            (T, n) == ("A", 1)

        return render_template("type.html", D=Dynkin(T, n))
    except ValueError:
        pass


@app.route("/<string:D>/<int:k>")
def show_grassmannian(D, k):
    T = D[0]
    n = int(D[1:])

    if T in ["B", "C"] and n == 1:
        return redirect("/A1/1")

    return render_template("grassmannian.html", G=Grassmannian(T, n, k))


@app.route("/horospherical")
def show_horospherical_table():
    return render_template("horospherical.table.html", horospherical=horospherical_varieties)


@app.route("/horospherical/<X>")
def show_horospherical(X):
    return render_template("horospherical.html", X=horospherical_varieties[X])

# to make the current year accessible in templates
@app.context_processor
def inject_now():
    return {"now": datetime.utcnow()}


"""
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
"""
