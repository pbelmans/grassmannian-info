from flask import Flask
from flask import render_template

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
        elif n == 2 and k == 2: return "\\mathbb{{P}}^{{3}}"
        else: return "\\OGr({},{})".format(k, 2*n + 1)

    elif T == "C":
        if k == 1: return "\\mathbb{{P}}^{{{}}}".format(2*n - 1)
        elif k == n:
            if n == 2: return "\\mathrm{{Q}}^{{3}}"
            else: return "\\LGr({},{})".format(k, 2*n)
        else: return "\\SGr({},{})".format(k, 2*n)

    elif T == "D":
        if k == 1: return "\\mathrm{{Q}}^{{{}}}".format(2*n - 2)
        else: return "\\OGr({},{})".format(k, 2*n)

    elif T == "E":
        if n == 6 and k == 1: return "\\mathbb{OP}^2"
        if n == 6 and k == 6: return "\\mathbb{OP}^{2,\\vee}"

    elif T == "G":
        if k == 1: return latex("B", 3, 1)

    return "\\mathrm{{{}}}{}/\\mathrm{{P}}{}".format(T, n, k)


# load the Sage-generated JSON file
with open("grassmannians.json") as f:
    data = json.load(f)

    for G in data:
        T = G["type"][0]
        n = int(G["type"][1:])
        k = G["parabolic"]

        # cutoff for now
        if T in "ABCD" and n > 7: continue

        G["latex"] = latex(T, n, G["parabolic"])

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

        if n not in grassmannians[T]:
            grassmannians[T][n] = {}

        grassmannians[T][n][k] = G



print(grassmannians)


@app.route("/")
def hello_world():
    return render_template("table.html", grassmannians=grassmannians)