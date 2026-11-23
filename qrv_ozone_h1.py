"""Exact Q(1) ozone theorem and V(1) upper bound.

Source: Gomez--Venegas supplied draft, tab:Mt, thm:viaA. These
relations are reconstructed from the normal-element localization,
not independently transcribed from Zhang--Zhang.
"""
from collections import defaultdict
import sympy as sp
from qrv_candidates import datum,candidates
from qrv_center_low import center,reducer,NAMES

def polynomial_image(poly,X,Y,f):
    maps={}
    for k,g in enumerate(NAMES):
        A=X if k<2 else Y
        i=k if k<2 else k-2
        bas=NAMES[:2] if k<2 else NAMES[2:]
        maps[g]={bas[j]:A[i,j] for j in range(2) if A[i,j]}
    norm=reducer(f);out=defaultdict(lambda:sp.Integer(0))
    for word,coef in poly.items():
        partial={():coef}
        for letter in word:
            new=defaultdict(lambda:sp.Integer(0))
            for w,c in partial.items():
                for g,v in maps[letter].items():new[w+(g,)]+=c*v
            partial=new
        for w,c in partial.items():
            for nw,nc in norm(w).items():out[nw]+=c*nc
    return {w:c for w,c in out.items() if c}

def vector_poly(bs,vs):return {w:vs[i] for i,w in enumerate(bs) if vs[i]}

def fixes(poly,X,Y,f):return polynomial_image(poly,X,Y,f)==poly

def Q_certificate():
    # Derived exact central element in B_{(1,2)} from a rational nullspace.
    bs,ker=center('Q',1,2)
    assert len(ker)==1
    z=vector_poly(bs,ker[0])
    assert z=={('x2','y2','y2'):-1,('x2','y1','y1'):1,
               ('x1','y2','y2'):1,('x1','y1','y2'):-2,('x1','y1','y1'):1}
    # The central subalgebra gives the necessary signed monomial candidates;
    # full central generation is NOT assumed.
    cc=candidates('Q')
    assert len(cc)==8
    filt=[(X,Y) for X,Y in cc if fixes(z,X,Y,'Q')]
    assert len(filt)==4
    # Conjugation by the normal x1 is an ozone automorphism.
    # Y x1 = x1 D Y, with D = [[0,1],[-1,0]].
    Gx=sp.diag(1,-1)
    Gy=sp.Matrix([[0,-1],[1,0]])
    assert any(X==Gx and Y==Gy for X,Y in filt)
    assert all(not r for r in __import__('general_quadratic').verify_general(datum('Q'),Gx,Gy))
    pow_set=[(Gx**k,Gy**k) for k in range(4)]
    assert all(any(X==A and Y==B for A,B in filt) for X,Y in pow_set)
    assert Gy**2==-sp.eye(2) and Gy**4==sp.eye(2)
    # Verify all known central constraints explicitly.
    tests=[(2,0),(0,2),(0,4),(4,0),(1,2)]
    for a,b in tests:
        bwords,space=center('Q',a,b)
        for v in space:
            poly=vector_poly(bwords,v)
            assert all(fixes(poly,X,Y,'Q') for X,Y in filt)
    return {'family':'Q','h':'1','upper_bound_before_new_center':len(cc),
            'upper_bound_after_new_center':len(filt),'ozone_group':'C4',
            'proof':'upper bound via exact central constraints and degree-one matrices; lower bound via Ad_x1',
            'new_central_element': '+'.join(str(c)+'*'+''.join(w) for w,c in z.items())}

def V_certificate():
    cc=candidates('V');assert len(cc)==4
    tests=[(1,0),(2,0),(0,2),(0,4),(2,2),(1,2),(3,2)]
    for a,b in tests:
        bwords,space=center('V',a,b)
        for v in space:
            poly=vector_poly(bwords,v)
            assert all(fixes(poly,X,Y,'V') for X,Y in cc)
    return {'family':'V','h':'1','upper_bound':len(cc),
            'status':'Only upper bound; the four candidates have been checked against selected central pieces, not the entire center',
            'tested_bidegrees':tests}

if __name__=='__main__':
 import json
 print(json.dumps({'Q':Q_certificate(),'V':V_certificate()},indent=2))
