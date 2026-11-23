"""Exact symbolic checks for quotient-conic correspondences of Q and V.

These checks prove polynomial identities for the dynamical reduction in the
Gomez--Venegas working manuscript. They DO NOT compute the ozone groups.
The R family needs a separate degree-two normal-element presentation.
"""
import sympy as sp

t, mu = sp.symbols('t mu')

def data(family):
    if family == 'Q':
        M = sp.Matrix([[-1, 1-t],[1+t, 1]])
        conic = mu**2 + t**2 - 2
        return M,conic,sp.Matrix([mu,-t]),sp.Matrix([0,0])
    if family == 'V':
        M = sp.Matrix([[1,t],[1-t,0]])
        conic = mu**2-mu+t**2-t
        return M,conic,sp.Matrix([mu,1-t]),sp.Matrix([sp.Rational(1,2)]*2)
    raise ValueError(f'No quotient-conic datum provided for {family}')

def iterate(phi, n):
    r = sp.Matrix([t,mu])
    for _ in range(n):
        r = sp.Matrix([sp.expand(v.subs({t:r[0],mu:r[1]}, simultaneous=True)) for v in phi])
    return r

def certificate(family):
    M,conic,phi,center=data(family)
    char=sp.factor(M.charpoly(mu).as_expr())
    assert sp.expand(char-conic)==0,(family,char,conic)
    updated=sp.expand(conic.subs({t:phi[0],mu:phi[1]}, simultaneous=True))
    assert sp.expand(updated-conic)==0,family
    assert iterate(phi,4)==sp.Matrix([t,mu]),family
    assert iterate(phi,2)!=sp.Matrix([t,mu]),family
    J=phi.jacobian([t,mu])
    assert J**2==-sp.eye(2) and J**4==sp.eye(2),family
    assert phi.subs({t:center[0],mu:center[1]},simultaneous=True)==center,family
    u,v=t-center[0],mu-center[1]
    quartic=sp.expand((u+sp.I*v)**4)
    # Under a quarter turn, the fourth power is invariant.
    assert sp.simplify(quartic.subs({t:phi[0],mu:phi[1]},simultaneous=True)-quartic)==0
    return {'family':family,'char_poly':str(char),'vieta_map':list(map(str,phi)),
            'order':4,'conic_invariant':True,'fixed_point':list(map(str,center)),
            'warning':'Not an ozone-group or center certificate'}

if __name__ == '__main__':
    import json
    print(json.dumps([certificate(f) for f in ('Q','V')],indent=2))
