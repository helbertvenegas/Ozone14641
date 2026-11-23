import unittest
from math import isqrt

ORDERS={'C*':12,'Q*':36,'R*':36,'V*':36}

def required_degree_factor(group_order):
    x=group_order; d=1; p=2
    while p*p<=x:
        if x%p==0:
            e=0
            while x%p==0: x//=p;e+=1
            d*=p**((e+1)//2)
        p+=1
    if x>1: d*=x
    return d

class Starred(unittest.TestCase):
    def test_required_factor(self):
        for name,order in ORDERS.items():
            self.assertEqual(required_degree_factor(order),6,name)
            for d in range(1,120):
                if d*d%order==0:self.assertEqual(d%6,0)

if __name__=='__main__':unittest.main()
