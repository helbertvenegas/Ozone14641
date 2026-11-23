"""Exact finite-group checks for D and G, Gomez--Venegas (2026).

Scope: G(p,f,1) for p of finite order n>=3, f!=0; D gives necessary
central constraints only. No numerical approximations or false exhaustive claims.
"""
from math import gcd
from collections import Counter
from sympy import symbols, factor, simplify

a,t,c,d,p,f,h = symbols('a t c d p f h', nonzero=True)

def mixed_obstructions(family):
    # In source generators (D) or shifted x2' generator (G), with
    # x1->a*x1, x2/x2'->t*x2/x2', y1->c*y1, y2->d*y2.
    if family not in ('D','G'):
        raise ValueError(family)
    return (factor(t*c-a*d), factor(t*d-a*c))

def product_orders(pair):
    """Order of element of Z/n x Z/n x Z/2 using modular addition."""
    n, ka, kc, eps = pair
    m=1
    while (ka*m)%n or (kc*m)%n or (eps*m)%2:
        m += 1
    return m

def g_group(n):
    """Enumerate ozone candidates: x1->a*x1, x2'->(-1)^eps*a*x2',
       y1->c*y1, y2->(-1)^eps*c*y2.
       Scalars a,c are nth roots, represented by exponents mod n.
    """
    if n < 3:
        raise ValueError('n>=3')
    k=n//2
    elems=[]
    for ka in range(n):
        for kc in range(n):
            for eps in range(2):
                if n%2:
                    ok=(eps==0)  # u_n fixes sign
                else:
                    # ((-1)^eps a)^(n/2) = 1
                    # a=exp(2pi i ka/n): a^(n/2)=(-1)^ka
                    ok=((ka+eps*k)%2==0)
                if ok:
                    elems.append((n,ka,kc,eps))
    return elems

def order_spectrum(n):
    return dict(sorted(Counter(product_orders(e) for e in g_group(n)).items()))

def expected_factors(n):
    if n%2 or (n//2)%2:
        return (n,n)
    return (2,n//2,n)

def factor_spectrum(factors):
    from itertools import product
    from math import lcm
    hist=Counter()
    for exps in product(*[range(m) for m in factors]):
        ords=[m//gcd(m,k) for m,k in zip(factors,exps)]
        order=1
        for z in ords: order=lcm(order,z)
        hist[order]+=1
    return dict(sorted(hist.items()))

def d_necessary(n):
    """Count candidates constrained ONLY by four central elements in draft.
       This is an upper bound, NOT the full ozone group.
    """
    e=n if n%2==0 else 2*n
    count=0
    for ka in range(e):
        for kc in range(e):
            for eps in range(2):
                # c^n * (-1)^(eps*n)=1 as e-th roots
                if (kc*n + eps*n*(e//2)) % e == 0:
                    count+=1
    return count

def run():
    assert mixed_obstructions('D') == mixed_obstructions('G')
    for n in range(3,21):
        elems=g_group(n)
        assert len(elems)==n*n,(n,len(elems))
        assert order_spectrum(n)==factor_spectrum(expected_factors(n)),n
    # D(4,1) is a sharp warning: central subset admits 32 candidates,
    # but manuscript table reports only 16.
    assert d_necessary(4)==32
    assert d_necessary(3)==36
    print('PASS: G(h=1) exact order spectra and invariant factors, n=3..20.')
    print('G examples: n=3',expected_factors(3),'n=4',expected_factors(4),
          'n=6',expected_factors(6),'n=8',expected_factors(8))
    print('D necessary candidate counts: n=3',d_necessary(3),
          'n=4',d_necessary(4),'(NOT an exhaustive ozone count)')

if __name__=='__main__': run()
