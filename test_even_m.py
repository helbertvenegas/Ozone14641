import unittest
from two_parameter_even_m import invariant_factors, FAMILIES

class SixFamiliesEvenM(unittest.TestCase):
    def test_all_six_even_m(self):
        for f in FAMILIES:
            for m in (2,4,6,8,10):
                for mp in (1,2,3,4,5,8):
                    with self.subTest(f=f,m=m,mp=mp):
                        item=invariant_factors(m,mp,f)
                        self.assertEqual(item['invariant_factors'], (2,2*mp,2*mp))
                        self.assertEqual(item['order'],8*mp*mp)
    def test_historical_example(self):
        for f in FAMILIES:
            self.assertEqual(invariant_factors(2,2,f)['invariant_factors'],(2,4,4))
    def test_scope_guard(self):
        with self.assertRaises(ValueError): invariant_factors(3,2,'O')
        with self.assertRaises(ValueError): invariant_factors(2,2,'G')

if __name__=='__main__': unittest.main()
