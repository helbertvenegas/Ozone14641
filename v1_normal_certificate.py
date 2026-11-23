"""Exact certificate for Oz(V(h=1)) = C4.

All operations are exact SymPy rational arithmetic. This proves the claim
relative to the defining presentation reconstructed from the supplied draft's
normal-localization table (tab:Mt). It is not a generic-h theorem.
"""
import json
import sympy as sp
from qrv_candidates import datum, candidates
from qrv_center_low import normal_intertwiners, reducer
from general_quadratic import verify_general


def certificate():
    X = sp.Matrix([[-1, 1], [0, 1]])
    Y = sp.Matrix([[0, -1], [1, 0]])
    assert all(not r for r in verify_general(datum('V'), X, Y)), 'Not an automorphism of the presented algebra'
    basis, nullspace = normal_intertwiners('V', 1, 2, X, Y)
    assert len(nullspace) == 1
    w = {word: nullspace[0][i] for i, word in enumerate(basis) if nullspace[0][i]}
    expected = {
        ('x1', 'y1', 'y1'): sp.Integer(1),
        ('x1', 'y2', 'y2'): sp.Integer(-1),
        ('x2', 'y1', 'y1'): sp.Integer(-1),
        ('x2', 'y1', 'y2'): sp.Integer(-1),
    }
    assert w == expected
    # Independent of the nullspace routine, check all 4 intertwining identities
    norm = reducer('V')
    gen = ('x1', 'x2', 'y1', 'y2')
    matrices = [(X, 0), (X, 1), (Y, 0), (Y, 1)]
    for g, (M, r) in zip(gen, matrices):
        left = {}
        right = {}
        for word, coeff in w.items():
            for mon, val in norm(word + (g,)).items():
                left[mon] = left.get(mon, 0) + coeff*val
            for c in range(2):
                if M[r, c] == 0:
                    continue
                new_g = gen[c] if g[0] == 'x' else gen[2+c]
                for mon, val in norm((new_g,) + word).items():
                    right[mon] = right.get(mon, 0) + coeff*M[r, c]*val
        assert {mon: v for mon, v in left.items() if v} == {mon: v for mon, v in right.items() if v}, g
    assert X**2 == sp.eye(2)
    assert Y**2 == -sp.eye(2)
    assert Y**4 == sp.eye(2)
    # Upper bound is the 4 candidates coming from the known central elements
    # and preservation of all quadratic relations. It is not a center-generation claim.
    admissible = candidates('V')
    assert len(admissible) == 4
    assert all(any(A == X**j and B == Y**j for A,B in admissible) for j in range(4))
    assert len({tuple(X**j) + tuple(Y**j) for j in range(4)}) == 4
    return {
        'family': 'V', 'h': '1', 'result': 'Oz(V(1)) = C4',
        'method': 'upper bound via necessary central constraints; lower bound via normal-inner generator',
        'normal_element': 'x1*(y1^2-y2^2)-x2*(y1^2+y1*y2)',
        'normal_bidegree': [1,2], 'intertwining_identities_checked': 4,
        'quadratic_relation_residuals': 0, 'upper_bound': 4,
        'distinct_powers_of_normal_conjugation': 4,
        'source_limitation': 'defining relations reconstructed from supplied manuscript tab:Mt; generic twists not proved',
    }

if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, ensure_ascii=False))
