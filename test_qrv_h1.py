"""Focused regression tests for three h=1 ozone-group proofs."""
import unittest
from qrv_ozone_h1 import Q_certificate
from r1_normal_certificate import certificate as R_certificate
from v1_normal_certificate import certificate as V_certificate

class ThreeOzoneGroupsAtH1(unittest.TestCase):
    def test_Q1_cyclic_four(self):
        r = Q_certificate()
        self.assertEqual(r['ozone_group'], 'C4')
        self.assertEqual(r['upper_bound_after_new_center'], 4)

    def test_R1_cyclic_four(self):
        r = R_certificate()
        self.assertEqual(r['candidates_after_new_center'], 4)
        self.assertEqual(r['distinct_powers_of_normal_conjugation'], 4)

    def test_V1_cyclic_four(self):
        r = V_certificate()
        self.assertEqual(r['upper_bound'], 4)
        self.assertEqual(r['distinct_powers_of_normal_conjugation'], 4)
        self.assertEqual(r['intertwining_identities_checked'], 4)

if __name__ == '__main__':
    unittest.main()
