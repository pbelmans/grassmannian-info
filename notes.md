# Ideas

* Dynkin diagrams: switch between Dynkin diagram and name? display both? only on hover?

* explain the coincidences: visual indication
  - accidental isomorphisms: link to the G/P for which Aut(G/P)=G
  - Dynkin diagram symmetries
  the visual indication should be an exclamation mark?

* on the dedicated page: give the full name such as
  - Lagrangian Grassmannian
  - Freudenthal variety
  - Cayley plane

* visual indication of type of collection:
  - use a 2x3 box, with borders
  - first row: conjecture, second row: proof
  - first column: exceptional collection, second column: Lefschetz collection, third column: minimal Lefschetz collection
  - encoded via six bits

* use "good" names
  - P^n
  - Q^n
  - LGr(n,2n) instead of SGr(n,2n)
  - ...?

* have a dedicated name function that deals with all the special cases (!)

* overview of references: table of the same format as on the overview table:

            | collection | Lefschetz | minimal Lefschetz
  ------------------------------------------------------
  candidate |
      proof |

  - make it zebra, with candidate in gray
  - in each cell: list of references, with MR or arXiv identifier
  - below the table: list of the references
  - hovering in the table highlights the reference below and vice versa


# To be done

* create a file with the numerical data
* create a file which takes

# Data format for the numerical data
Use JSON output from Python in Sage for this.

```json
[
  {
    "D" : "A4",
    "k" : 1,

    "dimension" : 4,
    "index" : 5,
    "K_0" : 5
  },
]
```

* Use this to initialise a  `grassmannians` variable with all cases under consideration.
* Make this into a `Grassmannian` class:
  - implement naming in Python

  - implement "alternatives" in Python

  - implement the bitmap in Python
    = set the following fields to an empty list []

      Grassmannian.collection.candidate
      Grassmannian.collection.proof
      Grassmannian.lefschetz.candidate
      Grassmannian.lefschetz.proof
      Grassmannian.minimal.candidate
      Grassmannian.minimal.proof

  - implement the logic in Jinja: if G.collection.candidate or G.collection.proof


# Data format for the exceptional collections

Via a Python file: use a template like

```python
reference = "MR"
reference = "arXiv"
for k in range(1, 11):
  grassmannians["A"][k].collection.candidate.append(reference)
```

Eventually one could turn this into a more elaborate system, with additional comments etc.
