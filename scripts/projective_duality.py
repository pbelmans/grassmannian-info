"""Compute projective-dual degrees from Tevelev, Theorem 7.1.

Run with Sage's Python: sage -python scripts/projective_duality.py
The ample generator is omega_k, including the spinor embeddings.
All arithmetic in the localization sum is integral; no numerical rounding.
"""

import argparse
import json
from math import prod
from pathlib import Path
from time import perf_counter

import numpy as np
from sage.all import RootSystem, ZZ, lcm, matrix, vector


ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / 'references/projective-duality.json'


def weight_orbit(simple_roots, k):
    """W/W_P as the orbit of omega_k in fundamental-weight coordinates."""
    highest = tuple(int(i == k - 1) for i in range(len(simple_roots)))
    weights, seen = [highest], {highest}
    for weight in weights:
        for i, root in enumerate(simple_roots):
            if not weight[i]:
                continue
            reflected = tuple(a - weight[i] * b for a, b in zip(weight, root))
            if reflected not in seen:
                seen.add(reflected)
                weights.append(reflected)
    return weights


def geometry(T, n, k):
    """Dual defect in the ample-generator embedding: Tevelev, Theorem 9.19."""
    empty = ((T == 'A' and k in {1, n}) or (T == 'C' and k == 1)
             or (T, n, k) == ('B', 2, 2)
             or (T == 'D' and n == 3 and k in {2, 3}))
    defect = 0
    if T == 'A' and n >= 4 and n % 2 == 0 and k in {2, n - 1}:
        defect = 2
    if (T, n, k) in {('B', 4, 4), ('D', 5, 4), ('D', 5, 5)}:
        defect = 4
    return empty, defect


def localization(roots, coroots, simple_roots, weights, dimension, defect, parameters):
    """Integrate c_(d-delta)(J_1(O(1))) H^delta by fixed-point localization.

    At the point indexed by mu=w(omega_k), H has weight -mu. Cotangent
    weights are the roots alpha with <mu, alpha^vee> > 0. The jet weights
    are -mu and alpha-mu. Thus c_(d-j)(J_1) is the coefficient of t^(j+1)
    in (t-mu) product_alpha(t+alpha-mu).

    Evaluate characters on a regular integral cocharacter. The product of
    all positive-root values clears every Euler-class denominator, so the
    entire calculation uses integer additions, multiplications and exact
    divisions. Also integrate H^d to independently check the polarization.
    """
    point = matrix(ZZ, simple_roots).solve_right(vector(ZZ, parameters))
    point *= lcm(x.denominator() for x in point)
    point = tuple(int(x) for x in point)
    simple_values = [sum(a * b for a, b in zip(root, point)) for root in simple_roots]
    values = [sum(a * b for a, b in zip(root, simple_values)) for root in roots]
    assert all(value > 0 for value in values)
    common = prod(values)
    totals = [0] * (defect + 1)
    degree_total = 0
    coroot_matrix = np.asarray(coroots, dtype=np.int64).T

    for start in range(0, len(weights), 2048):
        batch = weights[start:start + 2048]
        pairings = np.asarray(batch, dtype=np.int64) @ coroot_matrix
        for weight, pairing in zip(batch, pairings.tolist()):
            mu = sum(a * b for a, b in zip(weight, point))
            cotangent = [value if sign > 0 else -value
                         for value, sign in zip(values, pairing) if sign]
            assert len(cotangent) == dimension
            denominator = prod(-value for value in cotangent)
            quotient, remainder = divmod(common, denominator)
            assert remainder == 0
            degree_total += quotient * (-mu) ** dimension
            # Only the final delta+2 coefficients are needed, even for E8.
            coefficients = [-mu, 1] + [0] * defect
            for alpha in cotangent:
                value = alpha - mu
                for j in range(defect + 1, 0, -1):
                    coefficients[j] = coefficients[j] * value + coefficients[j - 1]
                coefficients[0] *= value
            for j in range(defect + 1):
                totals[j] += quotient * coefficients[j + 1] * (-mu) ** j

    integrals = []
    for total in [degree_total, *totals]:
        quotient, remainder = divmod(total, common)
        assert remainder == 0, 'localization did not give an integer'
        integrals.append(quotient)
    return integrals


def compute(row):
    D, k = row['type'], row['parabolic']
    T, n = D[0], int(D[1:])
    R = RootSystem([T, n])
    positive = list(R.root_lattice().positive_roots())
    roots = [tuple(int(x) for x in alpha.to_vector()) for alpha in positive]
    coroots = [tuple(int(x) for x in alpha.associated_coroot().to_vector()) for alpha in positive]
    simple_roots = [tuple(int(x) for x in alpha.to_vector())
                    for alpha in R.weight_lattice().simple_roots().values()]
    dimension = sum(root[k - 1] != 0 for root in roots)
    weights = weight_orbit(simple_roots, k)
    empty, defect = geometry(T, n, k)
    first = localization(roots, coroots, simple_roots, weights, dimension, defect, [1] * n)
    second = localization(roots, coroots, simple_roots, weights, dimension, defect, list(range(1, n + 1)))
    assert first == second, (D, k, 'cocharacter dependence')
    assert first[0] == row['degree'], (D, k, 'incorrect degree of the Grassmannian')
    assert all(value == 0 for value in first[1:-1]), (D, k, 'incorrect dual defect')
    assert first[-1] == 0 if empty else first[-1] > 0, (D, k, 'incorrect empty dual')
    return {
        'dimension': dimension, 'fixed_points': len(weights),
        'degree_of_grassmannian': first[0],
        'codimension': None if empty else defect + 1,
        'degree': None if empty else first[-1],
        'jet_integrals': first[1:],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='recompute and compare with the stored values')
    parser.add_argument('--case', help='compute a single case, e.g. E8/4, without changing the data file')
    args = parser.parse_args()
    data = json.loads((ROOT / 'data/grassmannians.json').read_text())
    rows = [row for row in data if int(row['type'][1:]) <= 8]
    results = {}
    for row in rows:
        key = f"{row['type']}/{row['parabolic']}"
        if args.case and key != args.case:
            continue
        started = perf_counter()
        results[key] = compute(row)
        print(key, results[key], f'({perf_counter() - started:.2f}s)', flush=True)
    if args.case:
        assert results, f'No catalogue entry: {args.case}'
        return
    if args.check:
        assert results == json.loads(TARGET.read_text())
        print(f'All {len(results)} stored records agree with both exact evaluations.')
    else:
        TARGET.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
