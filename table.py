from flask import Flask
from flask import render_template, redirect

import json

app = Flask(__name__)

grassmannians = {letter : {} for letter in "ABCDEFG"}

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
        elif n == 4 and k in [3, 4]: return latex("D", 4, 1)
        elif k <= n - 2: return "\\OGr({},{})".format(k, 2*n)
        else: return "\\OGr_+({},{})".format(n, 2*n)

    elif T == "E":
        if n == 6 and k == 1: return "\\mathbb{OP}^2"
        if n == 6 and k == 6: return "\\mathbb{OP}^{2,\\vee}"

    elif T == "G":
        if k == 1: return latex("B", 3, 1)

    return "\\mathrm{{{}}}_{{{}}}/\\mathrm{{P}}_{{{}}}".format(T, n, k)


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


