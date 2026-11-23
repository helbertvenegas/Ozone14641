"""Exact degree-two relation verifier for Zhang--Zhang monomial presentations.
This checks algebra-endomorphism relations only; it does not prove center fixing or ozone completeness.
Authors of underlying mathematics: James Gómez and Helbert Venegas.
"""
import sympy as S
from itertools import product
I=S.I

def family_data(name):
    """Concrete PI specialization of each of 14 monomial families, from tab:mono.
    Coordinate plane parameter p or q is interpreted in the chosen monomial basis.
    """
    h=(1+I)/2
    f=(1+I)/2; g=(1-I)/2
    data={
      'B':('I',('coord','coord'),(I,I),[1,1,-1,1],I),
      'E':('I',('balanced','coord'),(None,I),[1+I,1-I,-1+I,-1-I],S.sqrt(-I/2)),
      'J':('I',('coord','balanced'),(I,None),[1+I,-1+I,1-I,-1-I],S.sqrt(S.Rational(1,2))),
      'K':('0',('coord','balanced'),(1,None),[1,I,1,-I],1),
      'L':('II',('coord','coord'),(1,-1),[-1,1,-1,1],1),
      'M':('I',('balanced','balanced'),(None,None),[I*(I+1),-I*(I-1),I*(I-1),-I*(I+1)],h),
      'N':('III',('coord','balanced'),(-1,None),[f-g,f+g,-f-g,g-f],1),
      'O':('I',('balanced','balanced'),(None,None),[1-I,1+I,1+I,1-I],h),
      'P':('II',('balanced','coord'),(None,-1),[1+I,1-I,1-I,1+I],h),
      'S':('0',('balanced','balanced'),(None,None),[2,-2,-2,-2],S.Rational(1,2)),
      'T':('II',('balanced','balanced'),(None,None),[2*I,-2,-2*I,-2],S.Rational(1,2)),
      'U':('III',('balanced','balanced'),(None,None),[2*I,-2*I,-2,-2],S.Rational(1,2)),
      'W':('III',('coord','balanced'),(1,None),[I*(1-I),-I*(1+I),I*(1+I),-I*(1-I)],h),
      'Z':('II',('balanced','coord'),(None,1),[I*(1-I),-I*(1+I),-I*(1+I),I*(1-I)],h),
    }
    typ,sides,q,ls,tw=data[name]
    return {'family':name,'type':typ,'sides':sides,'q':q,'lambda':tuple(S.simplify(v*(1 if name=='N' else tw)) for v in ls),'twist':tw}

# Words are strings (e1,e2,z1,z2), mapping to free degree-two monomials.
E=('e1','e2'); Z=('z1','z2'); GEN=E+Z

def add(poly,word,coeff):
    poly[word]=S.simplify(poly.get(word,0)+coeff)

def multiply_linear(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():add(out,(x,y),c*d)
    return out

def reduced_word(word,d):
    a,b=word
    if a in Z and b in E:
        ia=Z.index(a);ib=E.index(b)
        typ=d['type'];js=1-ib if typ in ('I','III') else ib; jt=1-ia if typ in ('I','II') else ia
        return ((E[js],Z[jt]),d['lambda'][2*ia+ib])
    if a=='e2' and b=='e2' and d['sides'][0]=='balanced':return (('e1','e1'),S.Integer(1))
    if a=='z2' and b=='z2' and d['sides'][1]=='balanced':return (('z1','z1'),S.Integer(1))
    if a=='e2' and b=='e1' and d['sides'][0]=='coord':return (('e1','e2'),d['q'][0])
    if a=='z2' and b=='z1' and d['sides'][1]=='coord':return (('z1','z2'),d['q'][1])
    return (word,S.Integer(1))

def reduce_poly(poly,d):
    out={}
    for word,c in poly.items():
        w,k=reduced_word(word,d);add(out,w,c*k)
    return {k:S.simplify(v) for k,v in out.items() if S.simplify(v)!=0}

def relation_pairs(d):
    pairs=[]
    for block,side,q in zip((E,Z),d['sides'],d['q']):
        if side=='balanced':pairs.append(({(block[1],block[1]):1},{(block[0],block[0]):1}))
        else:pairs.append(({(block[1],block[0]):1},{(block[0],block[1]):q}))
    for a in range(2):
        for b in range(2):
            w=(Z[a],E[b]);rw,k=reduced_word(w,d)
            pairs.append(({w:1},{rw:k}))
    return pairs

def transform_poly(p,images):
    out={}
    for (a,b),k in p.items():
        for word,v in multiply_linear(images[a],images[b]).items():add(out,word,k*v)
    return out

def verify_linear_action(d,X,Y):
    """Returns six exact reduced degree-two residuals, all must vanish.
    Matrix convention: row i contains coefficients of image(generator i).
    """
    matrices=(S.Matrix(X),S.Matrix(Y))
    if any(m.shape!=(2,2) for m in matrices):raise ValueError('each block must be 2x2')
    images={}
    for block,m in zip((E,Z),matrices):
        for i in range(2):images[block[i]]={block[j]:m[i,j] for j in range(2)}
    residuals=[]
    for lhs,rhs in relation_pairs(d):
        p=transform_poly(lhs,images);q=transform_poly(rhs,images)
        for k,v in q.items():add(p,k,-v)
        residuals.append(reduce_poly(p,d))
    return residuals

def invertible(X,Y):return S.simplify(S.det(S.Matrix(X))*S.det(S.Matrix(Y)))!=0

def diagonal_sign_solutions(d):
    sols=[]
    for signs in product((-1,1),repeat=4):
        X=S.diag(*signs[:2]);Y=S.diag(*signs[2:]);res=verify_linear_action(d,X,Y)
        if all(not x for x in res):sols.append(signs)
    return sols

def permutation_candidates(d):
    ident=S.eye(2); swap=S.Matrix([[0,1],[1,0]])
    return {f'{a}_{b}':all(not p for p in verify_linear_action(d,X,Y)) for a,X in [('id',ident),('swap',swap)] for b,Y in [('id',ident),('swap',swap)]}

if __name__=='__main__':
    for fam in 'BEJKLM NOPSTUWZ'.replace(' ',''):
        d=family_data(fam)
        assert len(verify_linear_action(d,S.eye(2),S.eye(2)))==6
        assert all(not r for r in verify_linear_action(d,S.eye(2),S.eye(2)))
        ss=diagonal_sign_solutions(d)
        assert ss
        print(fam,'signed endomorphisms:',len(ss),'pure permutation patterns:',permutation_candidates(d))
