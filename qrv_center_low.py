"""Exact low-bidegree center computation in Q,V at h=1.
Not an exhaustive center presentation.
"""
from functools import lru_cache
from collections import defaultdict
import sympy as sp
from qrv_candidates import datum
NAMES=('x1','x2','y1','y2')

def reducer(f):
 d=datum(f)
 rules={('x2','x1'):{('x1','x2'):d.q12,('x1','x1'):d.q11},('y2','y1'):{('y1','y2'):d.p12,('y1','y1'):d.p11}}
 for iy in range(2):
  for ix in range(2):
   row=2*iy+ix
   rules[(NAMES[2+iy],NAMES[ix])]={('x1','y1'):d.sigma[row,0],('x2','y1'):d.sigma[row,1],('x1','y2'):d.sigma[row,2],('x2','y2'):d.sigma[row,3]}
 @lru_cache(None)
 def norm(word):
  for i in range(len(word)-1):
   if word[i:i+2] in rules:
    out=defaultdict(lambda:sp.Integer(0))
    for pair,c in rules[word[i:i+2]].items():
     if c==0:continue
     for w,v in norm(word[:i]+pair+word[i+2:]).items():out[w]+=c*v
    return {w:v for w,v in out.items() if v!=0}
  return {word:sp.Integer(1)}
 return norm

def words(a,b):
 return [tuple(['x1']*i+['x2']*(a-i)+['y1']*j+['y2']*(b-j)) for i in range(a+1) for j in range(b+1)]

def center(f,a,b):
 norm=reducer(f); bs=words(a,b); cols=[]; rowkeys=set()
 for z in bs:
  col={}
  for g in NAMES:
   left=norm(z+(g,));right=norm((g,)+z)
   for w in set(left)|set(right):
    c=left.get(w,0)-right.get(w,0)
    if c:
     col[(g,w)]=c;rowkeys.add((g,w))
  cols.append(col)
 rows=sorted(rowkeys)
 if not rows:return bs,[sp.eye(len(bs))[:,i] for i in range(len(bs))]
 M=sp.SparseMatrix(len(rows),len(bs),{(i,j):col[row] for j,col in enumerate(cols) for i,row in enumerate(rows) if row in col})
 return bs,M.nullspace()

if __name__=='__main__':
 for f in ('Q','V'):
  print('Family',f,flush=True)
  for a,b in [(1,0),(0,2),(1,2),(2,0),(2,2),(3,0),(3,2),(1,4)]:
   try:
    bs,ker=center(f,a,b)
    print((a,b),'dim',len(ker),flush=True)
    for v in ker[:2]:print('   ', ' + '.join(f'{c}*{"".join(bs[j])}' for j,c in enumerate(v) if c),flush=True)
   except Exception as ex:print('ERROR',ex)

def normal_intertwiners(f,a,b,X,Y):
 norm=reducer(f);bs=words(a,b); images={}
 for k,g in enumerate(NAMES):
  if k<2: images[g]={NAMES[j]:X[k,j] for j in range(2) if X[k,j]}
  else:images[g]={NAMES[2+j]:Y[k-2,j] for j in range(2) if Y[k-2,j]}
 cols=[];keys=set()
 for z in bs:
  col={}
  for g in NAMES:
   left=norm(z+(g,));right=defaultdict(lambda:sp.Integer(0))
   for h,c in images[g].items():
    for w,v in norm((h,)+z).items():right[w]+=c*v
   for w in set(left)|set(right):
    c=left.get(w,0)-right.get(w,0)
    if c:col[(g,w)]=c;keys.add((g,w))
  cols.append(col)
 rows=sorted(keys)
 M=sp.SparseMatrix(len(rows),len(bs),{(i,j):col[row] for j,col in enumerate(cols) for i,row in enumerate(rows) if row in col})
 return bs,M.nullspace()
