"""Exact certificates for the additional small-semisimple entries.

Run: sage scripts/quantum_semisimplicity.sage
The quantum Chevalley formula is Fulton–Woodward, On the quantum product of
Schubert classes, J. Algebraic Geom. 13 (2004), Theorem 10.1.
All computations take place over QQ at q=1, without numerical eigenvalues.
"""

import json
import sys
from collections import deque
from pathlib import Path


def characteristic_polynomial(T, n, k):
    roots = RootSystem([T, n]).root_lattice()
    ambient = RootSystem([T, n]).ambient_space()
    W = ambient.weyl_group()
    omega = ambient.fundamental_weights()[k]
    simples = W.simple_reflections()
    # w(omega_k) identifies the right coset w W_P. Breadth-first search gives
    # its minimal-length representative without enumerating the whole W.
    representatives = {omega: W.one()}
    queue = deque([omega])
    while queue:
        weight = queue.popleft()
        for s in simples.values():
            v = s * representatives[weight]
            target = v.action(omega)
            if target not in representatives:
                representatives[target] = v
                queue.append(target)
    weights = list(representatives)
    indices = {weight: i for i, weight in enumerate(weights)}
    nonparabolic = [a for a in roots.positive_roots() if a.coefficient(k)]
    two_rho = sum(a.to_ambient() for a in nonparabolic)
    reflection_data = []
    for a in nonparabolic:
        av = a.to_ambient()
        coroot = 2 * av / av.inner_product(av)
        reflection = W.one()
        for j in a.associated_reflection():
            reflection *= simples[j]
        reflection_data.append((reflection, ZZ(two_rho.inner_product(coroot)),
                                QQ(omega.inner_product(coroot))))
    M = matrix(QQ, len(weights))
    lengths = {weight: w.length() for weight, w in representatives.items()}
    for weight, w in representatives.items():
        for reflection, degree, coefficient in reflection_data:
            target = (w * reflection).action(omega)
            if lengths[target] in {lengths[weight] + 1, lengths[weight] + 1 - degree}:
                M[indices[target], indices[weight]] += coefficient
    return M.charpoly()


def certificate(T, n, k):
    polynomial = characteristic_polynomial(T, n, k)
    assert polynomial.gcd(polynomial.derivative()) == 1, (T, n, k)
    return {"dimension": polynomial.degree(), "coefficients": [str(c) for c in polynomial.list()],
            "gcd_with_derivative": "1"}


if __name__ == '__main__':
    target = Path(__file__).resolve().parent.parent / 'references/quantum-certificates.json'
    if '--check' in sys.argv:
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        # Independent normalisations: projective space and the doubled zero
        # eigenvalue of IGr(2,6), described by Chaput–Manivel–Perrin, §7.
        assert R(characteristic_polynomial('A', 4, 1)) == x**5 - 1
        igr = R(characteristic_polynomial('C', 3, 2))
        assert igr.valuation() == 2
        assert igr.gcd(igr.derivative()) == x
        for key, data in json.loads(target.read_text()).items():
            polynomial = R(data['coefficients'])
            assert polynomial.degree() == data['dimension'], key
            assert polynomial.gcd(polynomial.derivative()) == 1, key
        print('Known examples and all nine stored certificates verified over QQ.')
        sys.exit()
    cases = [('B', 6, 4), ('B', 7, 4), ('B', 7, 6), ('C', 5, 3),
             ('C', 6, 3), ('C', 6, 5), ('E', 6, 3), ('F', 4, 2)]
    results = {}
    for T, n, k in cases:
        key = f'{T}{n}/{k}'
        results[key] = certificate(T, n, k)
        print(key, results[key]['dimension'], 'square-free', flush=True)
    results['E6/5'] = results['E6/3']  # Dynkin diagram involution
    target.write_text(json.dumps(results, indent=2) + '\n')
