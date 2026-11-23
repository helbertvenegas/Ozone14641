from odd_m_analysis import group_data, diagonal_data
from fractions import Fraction as F

def predict(fam,U,V):
    if fam in 'MNP':
        k=group_data(fam,U,V)
        nonzero=[x[-1] for row in k['cosets'].values() for x in row]
        assert nonzero
        return 'B' if any(any(j==F(1,2) for j in s) for s in nonzero) else 'A'
    if U%2==0:
        c=diagonal_data(fam,U,V)
        assert c
        return 'B' if any(any(j==F(1,2) for j in s) for values in c.values() for s in values) else 'A'
    if V%2==1:return 'E'
    k=group_data(fam,U,V)
    vecs=set(x[-1] for rows in k['cosets'].values() for x in rows)
    assert len(vecs)>=2 and (F(0),F(0)) not in vecs
    return 'D'

def expected(fam,U,V):
    uodd=U%2==1
    vodd=V%2==1
    if fam in 'MP':
        if uodd: return 'A' if vodd else 'B'
        return 'A' if V%4==2 else 'B'
    if fam=='N':
        if uodd:return 'A' if V%4==2 else 'B'
        return 'A' if vodd else 'B'
    if not uodd:return 'A' if vodd else 'B'
    return 'E' if vodd else 'D'


def test_odd_m_cases():
    # U odd, or U=2*m with m odd; V any positive order
    count=0
    for fam in ['M','N','P','O','W','Z']:
        for U in [1,2,3,5,6,7,10,11,14,15]:
            for V in range(1,13):
                got=predict(fam,U,V)
                want=expected(fam,U,V)
                assert got==want,(fam,U,V,got,want)
                count+=1
    assert count==720
    return count

if __name__=='__main__':
    print('PASS:',test_odd_m_cases(),'exact cyclotomic-phase parameter instances')
