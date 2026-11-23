from fractions import Fraction as F
from itertools import product
from math import lcm

# action on r^m,rho^m,X^m',Y^m' as (sign exponent, root of unity phase)
# table from Gomez-Venegas manuscript, Table 'act'
raw={
'M': {'e':[(-1,'1'),(-1,'u^-1'),(1,'1'),(1,'v^-1')], 'z':[(-1,'-u'),(-1,'1'),(1,'-v'),(1,'1')], 'S':[(-1,'1'),(-1,'-1'),(1,'1'),(1,'1')], 'squares':(1,-1),'type':'I'},
'N': {'e':[(1,'-1'),(1,'-u'),(1,'-1'),(1,'v^-1')], 'z':[(-1,'u^-1'),(-1,'1'),(1,'-v'),(1,'1')], 'S':[(-1,'-1'),(-1,'1'),(1,'1'),(1,'1')], 'squares':(-1,1),'type':'III'},
'P': {'e':[(-1,'1'),(-1,'u^-1'),(1,'1'),(1,'-v^-1')], 'z':[(1,'u'),(1,'-1'),(1,'v'),(1,'-1')], 'S':[(-1,'1'),(-1,'1'),(1,'1'),(1,'-1')], 'squares':(1,1),'type':'II'},
'O': {'e':[(-1,'1'),(-1,'u'),(1,'1'),(1,'v^-1')], 'z':[(-1,'u^-1'),(-1,'1'),(1,'v'),(1,'1')], 'S':[(-1,'1'),(-1,'1'),(1,'1'),(1,'1')], 'squares':(1,1),'type':'I'},
'W': {'e':[(1,'1'),(1,'u'),(1,'1'),(1,'v^-1')], 'z':[(-1,'-u^-1'),(-1,'1'),(1,'v'),(1,'1')], 'S':[(-1,'-1'),(-1,'1'),(1,'-1'),(1,'1')], 'squares':(-1,1),'type':'III'},
'Z': {'e':[(-1,'1'),(-1,'u^-1'),(1,'1'),(1,'v^-1')], 'z':[(1,'u'),(1,'1'),(1,'v'),(1,'1')], 'S':[(-1,'1'),(-1,'1'),(1,'1'),(1,'1')], 'squares':(1,1),'type':'II'}
}

def parse(s,U,V):
    out=F(0)
    for term in s.replace('-u','+-u').replace('-v','+-v').split('+'):
        if not term: continue
        if term=='1':continue
        if term=='-1':out+=F(1,2)
        elif term=='u':out+=F(1,U)
        elif term=='u^-1':out-=F(1,U)
        elif term=='-u':out+=F(1,2)+F(1,U)
        elif term=='-u^-1':out+=F(1,2)-F(1,U)
        elif term=='v':out+=F(1,V)
        elif term=='v^-1':out-=F(1,V)
        elif term=='-v':out+=F(1,2)+F(1,V)
        elif term=='-v^-1':out+=F(1,2)-F(1,V)
        else: raise ValueError((s,term))
    return out%1

def powers(action,U,V,m,mp):
    return tuple((s, (parse(expr,U,V)*n)%1) for (s,expr),n in zip(action,[m,m,mp,mp]))

def compose(A,B):
    # A composed B (on generators). A(B(x))
    return tuple((a[0]*b[0],(b[1]+b[0]*a[1])%1) for a,b in zip(A,B))

def signature(d): return tuple((s,p%1) for s,p in d)

def group_data(family,U,V):
    assert U>0 and V>0
    m=U if U%2 else U//2
    mp=V if V%2 else V//2
    assert m%2==1
    data=raw[family]
    E=powers(data['e'],U,V,m,mp);Z=powers(data['z'],U,V,m,mp)
    G={'id':((1,F(0)),)*4,'e':E,'z':Z,'ez':compose(E,Z)}
    assert compose(E,Z)==compose(Z,E),(family,U,V,'noncommuting K action')
    S=powers(data['S'],U,V,m,mp)
    N=4*mp**2
    kconds=[]
    for eps in (1,-1):
        ephase=F(0) if eps==1 else F(1,2)
        tx= ephase if data['type'] in ['I','II'] else F(0)
        ty= ephase if data['type'] in ['I','III'] else F(0)
        # central images of r^m, rho^m are inversion with phases S * chi
        pr=(S[0][1]+m*tx)%1
        pq=(S[1][1]+m*ty)%1
        for ja,jb in product(range(4*mp),repeat=2):
            pa,pb=F(ja,4*mp),F(jb,4*mp)
            L=tuple([(S[0][0],pr),(S[1][0],pq),
                     (1,(S[2][1]+2*mp*pa)%1),
                     (1,(S[3][1]+2*mp*pb)%1)])
            matches=[name for name,g in G.items() if g==L]
            if not matches or 'id' in matches:continue
            # For |G|=4, equality on K implies actual ozone membership
            sqx=(F(0) if data['squares'][0]==1 else F(1,2))+tx+2*pa
            sqy=(F(0) if data['squares'][1]==1 else F(1,2))+ty+2*pb
            sqx%=1;sqy%=1
            parity=((sqx*mp)%1,(sqy*mp)%1)
            kconds.append((matches[0],eps,ja,jb,parity))
    cosets={k:[] for k in G}
    for item in kconds:cosets[item[0]].append(item)
    return {'family':family,'U':U,'V':V,'m':m,'mp':mp,'Gsize':len(set(G.values())),'cosets':cosets,'all':kconds}

if __name__=='__main__':
    for U in [3,6,5,10]:
      for V in [1,2,3,4,5,6]:
       for fam in raw:
        x=group_data(fam,U,V)
        if x['Gsize']!=4:continue
        # List coset splitting obstructions as square parity of representative
        vals={cos:sorted(set(item[-1] for item in items)) for cos,items in x['cosets'].items() if items}
        print(fam,U,V,'mprime',x['mp'],'cosets', vals)

def diagonal_data(family,U,V):
 m=U if U%2 else U//2
 mp=V if V%2 else V//2
 assert m%2==1
 data=raw[family]
 E=powers(data['e'],U,V,m,mp);Z=powers(data['z'],U,V,m,mp)
 G={'id':((1,F(0)),)*4,'e':E,'z':Z,'ez':compose(E,Z)}
 out={}
 for eps in (-1,):
  ep=F(1,2)
  tx=ep if data['type'] in ['I','II'] else F(0)
  ty=ep if data['type'] in ['I','III'] else F(0)
  for ja,jb in product(range(4*mp),repeat=2):
   pa,pb=F(ja,4*mp),F(jb,4*mp)
   L=((1,m*tx%1),(1,m*ty%1),(1,2*mp*pa%1),(1,2*mp*pb%1))
   matches=[name for name,g in G.items() if g==L]
   if matches:
    par=(2*pa*mp%1,2*pb*mp%1)
    for n in matches:out.setdefault(n,set()).add(par)
 return out
