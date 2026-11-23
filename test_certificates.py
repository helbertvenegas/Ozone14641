"""Execute: python -m unittest discover -v"""
import unittest
import sympy as sp
from source_monomial import FAMILIES,valid_specialization,check_rebased_planes,crosscheck_actions,original_datum_from_monomial
from general_quadratic import verify_general
from run_certificates import run
class AllFamilies(unittest.TestCase):
    def test_parameter_exclusions_and_nonzero_mixed_coefficients(self):
        for f in FAMILIES:
            with self.subTest(family=f):valid_specialization(f)
    def test_change_of_basis_plane_relations(self):
        for f in FAMILIES:
            with self.subTest(family=f):self.assertTrue(check_rebased_planes(f))
    def test_quadratic_action_crosscheck(self):
        for f in FAMILIES:
            with self.subTest(family=f):self.assertEqual(crosscheck_actions(f),20)
    def test_original_quadratic_presentations(self):
        for f in FAMILIES:
            with self.subTest(family=f):
                self.assertTrue(all(not t for t in verify_general(original_datum_from_monomial(f),sp.eye(2),sp.eye(2))))
    def test_full_catalog_and_json(self):
        d=run()
        self.assertEqual(d['family_total'],26)
        self.assertEqual(d['sample_rows'],73)
        self.assertEqual(d['quadratic_action_crosschecks'],280)
        self.assertTrue(all(d['families'][f]['center_fixed'] is None for f in FAMILIES))
if __name__=='__main__':unittest.main()
