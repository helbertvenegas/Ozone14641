"""Quadratic PBW relation verifier for arbitrary trimmed double-Ore data.

Accepts q12,q11,p12,p11 and a 4x4 coefficient matrix Sigma where rows are
(y1*x1, y1*x2, y2*x1, y2*x2) and columns (x1*y1,x2*y1,x1*y2,x2*y2).
With a *correct* source transcription, applies to every Zhang-Zhang 14641 family.
This module verifies only algebra relations, not central action/completeness.
"""
from dataclasses import dataclass
import sympy as sp
from degree2_certificates import multiply_linear,add
X=('x1','x2');Y=('y1','y2')
@dataclass
class Datum:
    q12: object
    q11: object
    p12: object
    p11: object
    sigma: object
    def __post_init__(self):
        self.sigma=sp.Matrix(self.sigma)
        assert self.sigma.shape==(4,4)

def normal_form_degree_two(word,d):
    a,b=word
    if (a,b)==('x2','x1'):return {('x1','x2'):d.q12,('x1','x1'):d.q11}
    if (a,b)==('y2','y1'):return {('y1','y2'):d.p12,('y1','y1'):d.p11}
    if a in Y and b in X:
        r=2*Y.index(a)+X.index(b)
        return {('x1','y1'):d.sigma[r,0],('x2','y1'):d.sigma[r,1],('x1','y2'):d.sigma[r,2],('x2','y2'):d.sigma[r,3]}
    return {word:1}

def reduce_degree_two(poly,d):
    out={}
    for word,c in poly.items():
        for k,v in normal_form_degree_two(word,d).items():add(out,k,c*v)
    return {k:sp.simplify(v) for k,v in out.items() if sp.simplify(v)!=0}

def image(poly,maps):
    out={}
    for (a,b),k in poly.items():
        for word,c in multiply_linear(maps[a],maps[b]).items():add(out,word,k*c)
    return out

def verify_general(d,A,B):
    A,B=sp.Matrix(A),sp.Matrix(B)
    assert A.shape==(2,2) and B.shape==(2,2)
    maps={}
    for basis,m in ((X,A),(Y,B)):
        for i in range(2):maps[basis[i]]={basis[j]:m[i,j] for j in range(2)}
    rels=[({('x2','x1'):1},normal_form_degree_two(('x2','x1'),d)),
          ({('y2','y1'):1},normal_form_degree_two(('y2','y1'),d))]
    for y in Y:
        for x in X:rels.append(({(y,x):1},normal_form_degree_two((y,x),d)))
    residuals=[]
    for lhs,rhs in rels:
        out=image(lhs,maps)
        for word,c in image(rhs,maps).items():add(out,word,-c)
        residuals.append(reduce_degree_two(out,d))
    return residuals

def from_monomial_coordinate(d):
    """Cross-check only when both planes use their coordinate bases."""
    if d['sides']!=('coord','coord'):raise ValueError('Requires both sides coordinate')
    sigma=sp.zeros(4)
    for ia in range(2):
        for ib in range(2):
            typ=d['type'];jx=1-ib if typ in ('I','III') else ib;jy=1-ia if typ in ('I','II') else ia
            row=2*ia+ib;col=2*jy+jx
            sigma[row,col]=d['lambda'][row]
    return Datum(d['q'][0],0,d['q'][1],0,sigma)

if __name__=='__main__':
    from degree2_certificates import family_data,verify_linear_action
    for fam in ('B','L'):
        datum=from_monomial_coordinate(family_data(fam))
        for A,B in ((sp.eye(2),sp.eye(2)),(sp.diag(-1,1),sp.diag(1,-1))):
            g=verify_general(datum,A,B)
            m=verify_linear_action(family_data(fam),A,B)
            assert (all(not r for r in g))==(all(not r for r in m)),fam
        print(fam,'generic Sigma PBW reduction agrees with independent monomial engine')
