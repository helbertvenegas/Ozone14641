"""Exact necessary-condition candidate search for ozone automorphisms Q,V at h=1."""
import itertools
import sympy as s
from general_quadratic import Datum,verify_general

def datum(f,h=1):
 if f=='Q':
  n,xp='x1','x2'; D=h*s.Matrix([[0,1],[-1,0]]); M=s.Matrix([[-1,1-s.Symbol('t')],[1+s.Symbol('t'),1]]);q=-1
 elif f=='V':
  n,xp='x2','x1';D=h*s.eye(2);M=s.Matrix([[1,s.Symbol('t')],[1-s.Symbol('t'),0]]);q=1
 elif f=='R':
  # Independent published presentation (x3=y1,x4=y2), ZZ09 family R.
  # Sigma is parameterized by h and does NOT use the normal-localization model.
  return Datum(-1,0,-1,0,h*s.Matrix([[1,1,1,0],[0,0,1,0],[0,1,0,0],[0,-1,-1,1]]))
 else:raise ValueError(f)
 t=s.Symbol('t');A=M.diff(t);C=M.subs(t,0)
 sig=s.zeros(4)
 for j,y in enumerate(('y1','y2')):
  for k,x in enumerate(('x1','x2')):
   row=2*j+k
   if x==n:
    for iy in range(2):sig[row,2*iy+('x1','x2').index(n)]=D[j,iy]
   else:
    for iy in range(2):
     sig[row,2*iy+('x1','x2').index(xp)]=(D*A)[j,iy]
     sig[row,2*iy+('x1','x2').index(n)]+=(D*C)[j,iy]
 return Datum(q,0,-1,0,sig)

def monomials():
 for a,b in itertools.product((-1,1),repeat=2):
  yield s.diag(a,b)
  yield s.Matrix([[0,a],[b,0]])

def candidates(f):
 d=datum(f)
 if f=='Q':
  xs=list(s.diag(a,b) for a,b in itertools.product((-1,1),repeat=2))
 else:
  xs=[s.eye(2),s.Matrix([[-1,1],[0,1]])]
 results=[]
 for X,Y in itertools.product(xs,monomials()):
  residual=verify_general(d,X,Y)
  if all(not z for z in residual):results.append((X,Y))
 return results

if __name__=='__main__':
 for f in ('Q','V'):
  print('\n',f,datum(f).sigma)
  cc=candidates(f)
  print('size',len(cc))
  for x,y in cc:print('X=',list(x),'Y=',list(y))
