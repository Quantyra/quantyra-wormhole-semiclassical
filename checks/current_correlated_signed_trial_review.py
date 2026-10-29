"""Author review of the joint reference rejection on the original domain."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
import json,sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_correlated_signed_trial_test as m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
ROOT=m.ROOT;digest=lambda p:sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=m.STEM+'-attempt3');args=ap.parse_args()
source=ROOT/'evidence'/(args.receipt+'.json');raw=json.loads(source.read_bytes());ctx.prec=768
assert raw['status']=='specified_complete_joint_reference_certificate_rejected'
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_correlated_signed_trial_test.py').read_text(encoding='utf-8')
protected=json.loads(m.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
assert digest(Path(sys.executable))==rt['python_sha256']
assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
r=raw['result'];c={k:dec(v) for k,v in r['constants'].items()}
old=m.reuse_reference_trial()
assert r['numerical_field']==old['numerical_field']
assert r['whole_cell_queries']==old['whole_cell_queries']
assert r['polynomial_trial']==old['polynomial_trial']
Q=lambda x:Fraction(int(x['numerator']),int(x['denominator']))
v={k:Q(x) for k,x in r['exact_nominal_parameters'].items()}
N=[Q(x) for x in r['polynomial_trial']['numerator']]
D=[Q(x) for x in r['polynomial_trial']['denominator']]
mul,add,poles=m.polymul,m.polyadd,m.polypoles
lin=[v['Cv'],v['delta']];quad=[v['d0'],v['d1'],Fraction(1)]
def pp(p,n):
 out=[Fraction(1)]
 for _ in range(n):out=mul(out,p)
 return out
Nfirst=[v['K0']*v['kappa']**2*x for x in mul(pp(lin,2),[v['Cv']*v['lambda_n'],1])]
Dfirst=mul(poles([v[k] for k in ['muT','muS','lambda_l','lambda_n']]),pp(quad,2))
cj=-3*v['alpha']*v['kappa']**3/(4*v['gamma']*v['A']*v['At'])
Nj=[cj*x for x in mul(pp(lin,3),poles([v['muY'],v['muS']+2*v['scalar_c']]))]
Dj=mul(poles([v[k] for k in ['twice_a','muS','muT']]),pp(quad,3))
identity=add(mul(N,mul(Dfirst,Dj)),[-x for x in mul(D,add(mul(Nfirst,Dj),mul(Nj,Dfirst)))])
assert all(x==0 for x in identity),'Entire-spectrum signed rational identity'
assert len(N)==8 and len(D)==12 and all(x>0 for x in D)
assert N[-1]!=0 and v['delta']>0

# The averaging estimate is for the same reference after all three filters.
# It is not assumed to estimate the actual quantum/Riesz discrepancy.
average=r['averaged_reference_comparison']
assert [Q(x) for x in average['numerator']]==Nj
assert [Q(x) for x in average['denominator']]==Dj
from math import comb
g=min(x/Fraction(comb(9,j)) for j,x in enumerate(Dj));M=sum(abs(x) for x in Nj)
assert Q(average['positive_Bernstein_guard'])==g and g>0
assert Q(average['numerator_coefficient_sum'])==M
p,co,bs=m.fr.load();m.lr.model_rows(p,co);co=m.old.setup(p,co)
XH=m.qball(v['At'])*(8*p['gamma']/15).sqrt()*m.qball(M/g)
assert XH.overlaps(dec(average['complete_H_profile_X_upper']))
avg=c['original_X_to_Et_comparison']*arb('1e-110')*XH/(2*p['gamma'])
assert avg<arb('4.034e-24') and c['averaged_reference_energy_error']<arb('4.034e-24')

# A WHOLE actual physical interval, not a point or sampled quadrature, gives
# the lower bound. The complete quantum form's 1.01 comparison is retained.
lowers=[]
assert len(r['whole_cell_queries'])==4
for query,width in zip(r['whole_cell_queries'],['1e-4','1e-6','1e-8','1e-10']):
 dx=dec(query['x_width']);assert dx.contains(arb(width))
 xi=dec(query['xi']);assert xi.is_finite()
 assert all(dec(x).is_finite() for x in query['whole_cell_values'].values())
 l=arb(max(arb(0),abs(xi).lower()))*(p['gamma']*dx).sqrt()/arb('1.01')
 l-=c['complete_representation_energy_error']+c['averaged_reference_energy_error']
 lowers.append(l)
assert lowers[1]>arb('3.371e34')
assert c['complete_representation_energy_error']<arb('4.951e10')
qt=arb('.338')/(1-m.cs('current-filtered-mixed-green-residual')['epsilon'])
vcap=dec(m.lr.result('current-signed-adjoint-trial-attempt')['known_complete_certificate_debits']['psi_retained_energy_upper'])
alpha=m.cs('current-filtered-mixed-green-residual')['alpha_original_radius_cost']
quartic=alpha*(arb('.338')/arb('1.338'))**2*qt*vcap
# Independent multiplication need not reproduce an up()-selected endpoint
# bit for bit. Both outward expressions must prove the stated conservative
# bound, which is then subtracted again in the whole-cell lower test.
assert quartic<arb('6.428e31') and c['quartic_complete_energy_upper']<arb('6.428e31')
assert lowers[1]-quartic>arb('3.364e34')
assert c['complete_trial_energy_lower']>arb('3.364e34')
assert c['original_bolt_trial_report_lower']>arb('1.6877e-8')
assert c['trial_report_over_whole_allowance']>arb('10.124')
assert c['minimum_center_change_in_energy']>arb('3.032e34')
assert c['decisive_remaining_j_comparison_threshold']>arb('2.996e34')
assert c['proposed_half_gap_j_error_contract']>arb('1.49e34')
pref=m.cs('current-filtered-mixed-green-residual')['pref_t']
assert (pref*c['complete_trial_energy_lower']).overlaps(c['original_bolt_trial_report_lower'])
assert len(r['numerical_field']['upper_half_nodes'])==384 and len(r['contour_bounds'])==3
assert r['scope']['prescribed_averaged_reference_error_included']
assert r['scope']['specified_reference_certificate_rejected']
assert r['full_actual_j_response'] is None and r['actual_Zr'] is None
for key in ['actual_j_Riesz_or_filtered_response_evaluated','actual_signed_Zr_bounded','actual_bolt_screen_rejected',
 'original_admission','nonlinear_existence','exact_compatibility','physical_matching','validity_noise',
 'all_sector_stability','construction_obstruction','goal_complete']:
 assert r['scope'][key] is False
history=[]
for pth in sorted((ROOT/'evidence').glob(m.STEM+'*.json')):
 if pth==source or 'verification' in pth.name:continue
 rr=json.loads(pth.read_bytes());history.append(dict(path=str(pth.relative_to(ROOT)),sha256=digest(pth),status=rr['status']))
out=ROOT/'evidence'/(args.receipt+'-verification.json');assert not out.exists()
record=dict(status='specified_joint_reference_estimator_failure_reviewed_actual_response_open',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(source),review_source_sha256=digest(Path(__file__)),
 frozen_inputs_checked=len(raw['input_sha256_after']),protected_paths_checked=len(protected),
 review_type='author ordinary domain argument, exact signed polynomial identity, whole-cell lower and complete budget',
 independent_specialist_review=False,lean_verification=False,exact_rational_identity_coefficients=len(identity),
 whole_physical_cells_retained=4,prescribed_averaging_error_charged=True,entire_quartic_term_retained=True,
 reference_certificate_rejected=True,actual_response_bounded=False,construction_obstruction=False,goal_complete=False,
 preserved_receipts=history)
out.write_bytes((json.dumps(record,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in record.items() if not isinstance(v,list)},indent=2))
