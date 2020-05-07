from flask import Flask
from flask import render_template
from flask_scss import Scss

import json

app = Flask(__name__)
Scss(app, static_dir="static", asset_dir="assets")
print("hi")

grassmannians = {letter : {} for letter in "ABCDEFG"}

def latex(T, n, k):
    if T == "A":
        if k == 1: return "\\mathbb{{P}}^{{{}}}".format(n)
        elif k == n: return "\\mathbb{{P}}^{{{},\\vee}}".format(n)
        else: return "\\Gr({},{})".format(k, n + 1)

    elif T == "B":
        if k == 1: return "\\mathrm{{Q}}^{{{}}}".format(2*n - 1)
        else: return "\\OGr({},{})".format(k, 2*n + 1)

    elif T == "C":
        if k == 1: return "\\mathbb{{P}}^{{{}}}".format(2*n - 1)
        elif k == n: return "\\LGr({},{})".format(k, 2*n)
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

        if n not in grassmannians[T]:
            grassmannians[T][n] = {}

        grassmannians[T][n][k] = G


print(grassmannians)


@app.route("/")
def hello_world():
    return render_template("table.html", grassmannians=grassmannians)
