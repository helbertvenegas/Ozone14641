"""Deterministic, exact (symbolic) regression report for type (14641).

Does NOT certify pointwise action on a generating set of Z(B), nor ozone completeness.
"""
import csv,json,os
from pathlib import Path
import sympy as sp
from degree2_certificates import family_data,verify_linear_action
from source_monomial import (FAMILIES,valid_specialization,check_rebased_planes,
    original_datum_from_monomial,crosscheck_actions)
from general_quadratic import verify_general
ROOT=Path(__file__).resolve().parent

def run():
    with (ROOT/'family_catalog.csv').open(encoding='utf8',newline='') as fh:
        cat=list(csv.DictReader(fh))
    with (ROOT/'ozone_table.csv').open(encoding='utf8',newline='') as fh:
        rows=list(csv.DictReader(fh))
    assert len(cat)==26 and len({r['family'] for r in cat})==26
    assert len(rows)==73
    families={}
    for f in FAMILIES:
        d=family_data(f)
        params=valid_specialization(f)
        assert check_rebased_planes(f)
        datum=original_datum_from_monomial(f)
        assert all(not r for r in verify_general(datum,sp.eye(2),sp.eye(2)))
        cross=crosscheck_actions(f)
        families[f]={
            'test_specialization':{k:str(v) for k,v in params.items()},
            'mixed_scalars_nonzero':True,
            'original_plane_relations':True,
            'generic_PBW_identity':True,
            'symbolic_action_crosschecks':cross,
            'source_basis':'GV monomial manuscript Table tab:mono',
            'original_sigma_provenance':'reconstructed from monomial presentation, not independently transcribed from ZZ09',
            'verified_scope':'quadratic relations and basis-change consistency',
            'center_fixed':None,
            'ozone_completeness':None
        }
    for row in cat:
        f=row['family']
        if f not in families:
            families[f]={
                'source_classification':row['PI_classification_status'],
                'degree_two_matrix_loaded':False,
                'center_fixed':None,
                'ozone_completeness':None
            }
    assert len(families)==26
    samples={}
    for r in rows:
        f=r['familia'];samples[f]=samples.get(f,0)+1
    result={
        'schema_version':3,
        'authors':['James Gómez','Helbert Venegas'],
        'family_total':26,
        'sample_rows':len(rows),
        'monomial_relation_families':len(FAMILIES),
        'quadratic_action_crosschecks':sum(families[f]['symbolic_action_crosschecks'] for f in FAMILIES),
        'sample_counts_by_family':samples,
        'families':families,
        'disclaimer':'Zero relation mismatches does not establish full ozone group. Center generators and exhaustive automorphism solutions remain absent.'
    }
    with (ROOT/'certificate_report.json').open('w',encoding='utf8') as fh:
        json.dump(result,fh,indent=2,ensure_ascii=False)
    return result

if __name__=='__main__':
    r=run()
    print('Families:',r['family_total'],'sample rows:',r['sample_rows'])
    print('Monomial families:',r['monomial_relation_families'])
    print('Exact quadratic-action crosschecks:',r['quadratic_action_crosschecks'])
    print('Center/completeness certificates: pending')
