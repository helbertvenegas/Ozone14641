"""Reproducible one-by-one audit, NOT a proof of generic ozone groups.
Source material: GV manuscripts (uploaded), and v7 project.
"""
import csv, json, re, traceback
from pathlib import Path
from collections import Counter
import sympy as sp
from degree2_certificates import family_data, verify_linear_action
from source_monomial import FAMILIES, valid_specialization, check_rebased_planes
from qrv_ozone_h1 import Q_certificate
from r1_normal_certificate import certificate as R_certificate
from v1_normal_certificate import certificate as V_certificate
from dg_exact import g_group, order_spectrum, factor_spectrum, expected_factors, d_necessary
ROOT=Path(__file__).resolve().parent
FULL=set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
NO_PI=set('AHXY')
MONO=set(FAMILIES)
H1={'Q':Q_certificate,'R':R_certificate,'V':V_certificate}

def cyclic_factors(s):
    """Parse the literal LaTeX patterns in the 73-row historical table only."""
    s=s.replace('\\mathbb Z/','Z/').replace('\\times','*').replace(' ','').replace('$','')
    s=s.replace('\\left','').replace('\\right','')
    out=[]
    pat=r'\(?\(?Z/(\d+)\)?\)?(?:\^\{?(\d+)\}?)?'
    for match in re.finditer(pat,s):
        m=int(match.group(1));n=int(match.group(2) or 1)
        out.extend([m]*n)
    return out

def audit_samples():
    with (ROOT/'ozone_table.csv').open(encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=73: raise ValueError(f'expected 73 source rows, found {len(rows)}')
    out=[]
    for idx,r in enumerate(rows,1):
        errors=[]
        factors=cyclic_factors(r['grupo'])
        if not factors: errors.append('group expression not parsed')
        product=sp.prod(factors)
        if int(product)!=int(r['orden']):errors.append('cyclic factors vs stated order')
        if r['rango'] not in ('---','—',''):
            rank=int(r['rango']);ratio=int(r['cociente'])
            if rank!=int(r['orden'])*ratio:errors.append('rank/order ratio incorrect')
            if not sp.integer_nthroot(rank,2)[1]:errors.append('rank is not a square')
            if rank%int(r['orden']):errors.append('ozone group order does not divide rank')
        out.append({'row':idx,'family':r['familia'],'parameters':r['parametros'],
            'twist':r['h'],'factors':factors,'order':int(r['orden']),
            'checks_passed':not errors,'errors':errors,
            'scope':'arithmetic consistency only; no full center or exhaustivity check'})
    return out

def audit_family(f):
    result={'family':f,'tests':{},'verified':'none','generic_ozone_proved_by_software':False}
    try:
        if f in NO_PI:
            result['verified']='PI obstruction attributed to manuscript, not independently re-proved'
            result['tests']['catalogued_as_non_PI']=True
        if f in MONO:
            d=family_data(f)
            valid_specialization(f)
            assert check_rebased_planes(f)
            assert all(not p for p in verify_linear_action(d,sp.eye(2),sp.eye(2)))
            result['tests'].update({'parameter_restrictions_sample':True,'change_of_basis_plane_relations':True,
                                    'six_quadratic_relations_identity':True})
            result['verified']='one exact representative; source-monomial relation consistency'
        if f in H1:
            r=H1[f]()
            result['tests']['h_equals_1_normal_or_center_certificate']=True
            result['tests']['h_equals_1_summary']=str(r)
            result['verified']='h=1 specialized ozone C4 under sourced defining relations'
        if f=='G':
            for n in (3,4,5,6,8,12):
                assert len(g_group(n))==n*n
                assert order_spectrum(n)==factor_spectrum(expected_factors(n))
            result['tests']['h_equals_1_conditional_group_structure']=True
            result['verified']='conditional h=1 general theorem; finite arithmetic checks n=3,4,5,6,8,12'
        if f=='D':
            assert d_necessary(4)==32
            result['tests']['necessary_candidates_n4']=32
            result['verified']='necessary conditions only, no complete group'
        result['passed']=True
    except Exception as exc:
        result['passed']=False
        result['error']=repr(exc)
    if f in set('CFI'):
        result['verified']='historical examples only; generic PI criterion and Oz unresolved'
    if f in set('QRV'):
        result['generic_twists']='not certified'
    if f in MONO:
        result['manuscript_general_order_formula']='theorem-dependent, not independently certified here'
    return result

def main():
    with (ROOT/'family_catalog.csv').open(encoding='utf8',newline='') as fh:
        cat=list(csv.DictReader(fh))
    assert len(cat)==26 and {c['family'] for c in cat}==FULL
    rows=audit_samples()
    families={f:audit_family(f) for f in sorted(FULL)}
    counts={f:sum(1 for r in rows if r['family']==f) for f in sorted(FULL)}
    report={'schema_version':8,'authors':['James Gómez','Helbert Venegas'],
      'basis':'four GV drafts and the v7 code, with no independent ZZ09 Sigma audit',
      'source_cases':26,'sample_rows':len(rows),'sample_arithmetic_passed':sum(r['checks_passed'] for r in rows),
      'family_tests_passed':sum(r['passed'] for r in families.values()),
      'families_with_pi_members':22,'families_with_monomial_general_order_theorem':14,
      'families_with_h1_specialized_C4_proof':3,
      'generic_ozone_certified_independently_for_all_PI_parameters':0,
      'monomial_families':sorted(MONO),'no_PI_families':sorted(NO_PI),
      'historical_sample_counts':counts,'families':families,'historical_rows':rows,
      'limitations':['14 monomial results depend on supplied manuscript theorems',
          'D, G general twists, Q R V general twists, C F I generic family are not complete',
          '4 starred cases lack demonstrated complete center/rank',
          'exact source Section 4 ZZ09 Sigma matrices have not all been independently transcribed',
          'passing a relation check or a numerical example does not establish all ozone automorphisms']}
    (ROOT/'AUDIT_26_FAMILIES.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    lines=['# Individual mathematical audit of all 26 families','','Test labels are deliberately narrow; **PASS does not mean the full ozone group is proved**.','','| Family | PI scope | Independently executed in this package | Full general ozone certificate |','|---|---|---|---|']
    for f,r in families.items():
        scope='none in working PI classification' if f in NO_PI else 'members with PI parameters'
        tests=', '.join(k for k,v in r['tests'].items() if v is True) or 'catalog only'
        full='Not applicable' if f in NO_PI else 'Not independently certified'
        lines.append(f'| {f} | {scope} | {tests} | {full} |')
    lines.extend(['','## Priorities after this audit','',
      '- Obtain original matrices Sigma from ZZ09 and compare every entry under the stated basis conventions.',
      '- For every PI parameter stratum, certify a full generating set of the center or a rigorous alternative upper bound.',
      '- Derive all automorphisms satisfying the six quadratic relations and central constraints, including non-diagonal solutions.',
      '- Resolve generic twists for D, G, Q, R, V; resolve C, F, I and four starred rows.',
      '- Add proof of generic invariant factors for M, N, O, P, W and Z; Theorem A supplies orders rather than all abstract structures.',
      '',f'**Historical sample row arithmetic:** {sum(r["checks_passed"] for r in rows)}/{len(rows)}.',
      f'**Family-level program runs:** {sum(r["passed"] for r in families.values())}/26 (some are catalog-only).',''])
    (ROOT/'AUDIT_26_FAMILIES.md').write_text('\n'.join(lines),encoding='utf8')
    print(f'26-family audit: {report["family_tests_passed"]}/26 passed, sample arithmetic {report["sample_arithmetic_passed"]}/73.')
    for f,r in families.items(): print(f, 'PASS' if r['passed'] else 'FAIL', '-',r['verified'])
    return report
if __name__=='__main__':main()
