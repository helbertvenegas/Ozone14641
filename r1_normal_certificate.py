"""Exact symbolic proof certificate for Oz(R(1)) = C4.

R relations are independently transcribed from the family R presentation at
https://4d-as.ncag.info/families/ore-r/ (accessed 2026-11-22).
Exact arithmetic in Q, conditional on correctness of that source presentation.
"""
import itertools
import json
import sympy as sp
from general_quadratic import verify_general
from qrv_candidates import datum, monomials
from qrv_center_low import center, normal_intertwiners, reducer
from qrv_ozone_h1 import fixes, vector_poly


def certificate():
    X = sp.Matrix([[0, 1], [-1, 0]])
    Y = sp.Matrix([[0, -1], [1, 0]])
    d = datum('R')
    assert all(not residual for residual in verify_general(d, X, Y))
    # Normal element of bidegree (1,1): u = x2*y1 + x1*y2.
    basis, nullspace = normal_intertwiners('R', 1, 1, X, Y)
    assert len(nullspace) == 1
    u = vector_poly(basis, nullspace[0])
    assert u == {('x2', 'y1'): 1, ('x1', 'y2'): 1}
    normal = reducer('R')
    names = ('x1', 'x2', 'y1', 'y2')
    for k, gen in enumerate(names):
        mat = X if k < 2 else Y
        idx = k if k < 2 else k - 2
        lhs = {}
        rhs = {}
        for word, c in u.items():
            for nw, v in normal(word + (gen,)).items():
                lhs[nw] = lhs.get(nw, 0) + c*v
            for j in range(2):
                if mat[idx,j] == 0:
                    continue
                image_gen = names[j] if k < 2 else names[2+j]
                for nw,v in normal((image_gen,) + word).items():
                    rhs[nw] = rhs.get(nw, 0) + c*mat[idx,j]*v
        assert {w:c for w,c in lhs.items() if c} == {w:c for w,c in rhs.items() if c}, gen
    assert X**2 == -sp.eye(2) and Y**2 == -sp.eye(2)
    assert X**4 == sp.eye(2) and Y**4 == sp.eye(2)

    # All ozone automorphisms preserve bigrading, hence they preserve
    # the x,y planes. Each -1 quantum plane has monomial graded
    # automorphisms. Pointwise fixing x1²+x2², y1²+y2² forces signs.
    rel_candidates = [(A,B) for A,B in itertools.product(monomials(),monomials())
                      if all(not r for r in verify_general(d,A,B))]
    assert len(rel_candidates) == 8
    b,ns = center('R',1,1)
    assert len(ns)==1
    z = vector_poly(b,ns[0])
    assert z == {('x2','y2'):-1,('x2','y1'):1,('x1','y2'):1,('x1','y1'):1}
    fixed = [(A,B) for A,B in rel_candidates if fixes(z,A,B,'R')]
    assert len(fixed)==4
    assert all(any(A==X**j and B==Y**j for A,B in fixed) for j in range(4))
    return {
        'family':'R', 'h':'1', 'result':'Oz(R(1)) = C4',
        'relation_candidates':8, 'candidates_after_new_center':4,
        'central_element':'x1*y1+x1*y2+x2*y1-x2*y2',
        'normal_element':'x2*y1+x1*y2',
        'intertwining_identities_checked':4,
        'distinct_powers_of_normal_conjugation':4,
        'status':'Exact h=1 theorem conditional on the published R presentation; no generic-h claim',
    }

if __name__=='__main__':
    print(json.dumps(certificate(), indent=2))
