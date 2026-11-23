"""Source-traceable monomial specializations from GV monomial manuscript,
Table tab:mono (18 Nov 2026). Exact arithmetic only.

Not a source for nonmonomial families and not a center/ozone certificate.
"""
import sympy as sp
from degree2_certificates import family_data
from general_quadratic import Datum,verify_general

I=sp.I
FAMILIES=tuple('BEJKLM NOPSTUWZ'.replace(' ',''))


def valid_specialization(name):
    """Parameters obey defining exclusions and all four mixed scalars nonzero."""
    d=family_data(name)
    assert all(sp.simplify(t)!=0 for t in d['lambda']), (name,'zero mixed scalar')
    params={
        'B':{'p':I,'h':I}, 'E':{'p':I,'h':sp.sqrt(-I/2)},
        'J':{'q':I,'h':sp.sqrt(sp.Rational(1,2))},
        'K':{'q':1,'f':-1,'h':1},'L':{'q':1,'f':-1,'h':1},
        'M':{'f':-1,'h':(1+I)/2},
        'N':{'f':(1+I)/2,'g':(1-I)/2},
        'O':{'f':-1,'h':(1+I)/2},'P':{'f':-1,'h':(1+I)/2},
        'S':{'h':sp.Rational(1,2)},'T':{'h':sp.Rational(1,2)},
        'U':{'h':sp.Rational(1,2)},
        'W':{'f':1,'h':(1+I)/2},'Z':{'f':1,'h':(1+I)/2}
    }[name]
    if name in ('M','O','P','W','Z'):
        assert sp.simplify(params['f'])!=0
    if name in ('M','O','P'): assert params['f']!=1
    if name=='W':assert params['f']!=-1
    if name=='Z':assert params['f'] not in (0,-1)
    if name=='N':assert sp.simplify(params['f']**2-params['g']**2)!=0
    return params


def original_basis_data(name):
    """Matrices E and Z: monomial generators as linear forms in original x,y.
    Parameter s is the selected square root f in Table tab:mono.
    """
    name=name.upper(); eye=sp.eye(2)
    plus=lambda t:sp.Matrix([[1,t],[1,-t]])
    p={
        'B':(eye,eye,I,I),
        'E':(plus(I),eye,-1,I),
        'J':(eye,plus(I),I,-1),
        'K':(eye,plus(-I),1,-1),  # s=i, s^-1=-i
        'L':(eye,eye,1,-1),
        'M':(plus(1),plus(I),-1,-1),
        'N':(eye,plus(1),-1,-1),
        'O':(plus(-1),plus(I),-1,-1),
        'P':(plus(I),eye,-1,-1),
        'S':(plus(1),plus(1),-1,-1),
        'T':(plus(1),plus(I),-1,-1),
        'U':(plus(I),plus(1),-1,-1),
        'W':(plus(I),plus(1),1,-1),
        'Z':(plus(1),plus(I),-1,1)
    }
    return p[name]


def original_datum_from_monomial(name):
    """Reconstruct Sigma in original coordinates through an exact basis change.
    Must compare to independent ZZ09 Sigma before considering it a certificate.
    """
    d=family_data(name)
    E,Z,q,p=original_basis_data(name)
    mono=sp.zeros(4)
    for a in range(2):
        for b in range(2):
            jx=(1-b) if d['type'] in ('I','III') else b
            jy=(1-a) if d['type'] in ('I','II') else a
            mono[2*a+b,2*jx+jy]=d['lambda'][2*a+b]
    sigma=sp.kronecker_product(Z,E).inv()*mono*sp.kronecker_product(E,Z)[:, [0,2,1,3]]
    sigma=sigma.applyfunc(sp.simplify)
    return Datum(q,0,p,0,sigma)


def check_rebased_planes(name):
    """Verify each monomial-basis plane relation inside original q-planes."""
    d=family_data(name)
    E,Z,q,p=original_basis_data(name)
    def product(A,i,j,qp):
        v0,v1=A[i,0],A[i,1];w0,w1=A[j,0],A[j,1]
        return sp.Matrix([v0*w0, v0*w1+qp*v1*w0,v1*w1])
    for A,side,qp,qmono in ((E,d['sides'][0],q,d['q'][0]),(Z,d['sides'][1],p,d['q'][1])):
        if side=='balanced':
            residual=product(A,1,1,qp)-product(A,0,0,qp)
        else:
            residual=product(A,1,0,qp)-qmono*product(A,0,1,qp)
        assert all(sp.simplify(v)==0 for v in residual),(name,side,residual)
    return True

if __name__=='__main__':
    from degree2_certificates import verify_linear_action
    for name in FAMILIES:
        params=valid_specialization(name)
        assert check_rebased_planes(name)
        datum=original_datum_from_monomial(name)
        assert all(not x for x in verify_general(datum,sp.eye(2),sp.eye(2)))
        print(name,'valid params',params,'nonzero Sigma entries',len([v for v in datum.sigma if v!=0]))

def crosscheck_actions(name):
    """All 16 sign combinations + 4 pure swaps: two independent PBW engines."""
    from itertools import product
    from degree2_certificates import verify_linear_action
    d=family_data(name)
    E,Z,_,_=original_basis_data(name)
    datum=original_datum_from_monomial(name)
    patterns=[sp.eye(2),sp.Matrix([[0,1],[1,0]])]
    tested=0
    candidates=[(sp.diag(*signs[:2]),sp.diag(*signs[2:])) for signs in product((-1,1),repeat=4)]
    candidates +=list(product(patterns,patterns))
    for A,B in candidates:
        m_ok=all(not x for x in verify_linear_action(d,A,B))
        AA=E.inv()*A*E
        BB=Z.inv()*B*Z
        g_ok=all(not x for x in verify_general(datum,AA,BB))
        assert m_ok==g_ok, (name,A,B,m_ok,g_ok)
        tested+=1
    return tested
