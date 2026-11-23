"""Fast tests for the 26-family audit (specializations ≠ generic theorems)."""
import unittest
from audit_all_families import main

class AuditAll26(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.report=main()
    def test_catalog(self):
        self.assertEqual(len(self.report['families']),26)
        self.assertEqual(self.report['families_with_pi_members'],22)
        self.assertEqual(len(self.report['no_PI_families']),4)
    def test_source_rows(self):
        self.assertEqual(self.report['sample_arithmetic_passed'],73)
    def test_family_status_truthful(self):
        self.assertFalse(any(r['generic_ozone_proved_by_software'] for r in self.report['families'].values()))
        self.assertTrue(all(r['passed'] for r in self.report['families'].values()))
    def test_h1_case_certificates(self):
        for f in 'QRV':
            self.assertTrue(self.report['families'][f]['tests']['h_equals_1_normal_or_center_certificate'])
if __name__=='__main__':unittest.main()
