"""Author review of both original residuals, composed errors and finite field."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
import json,sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_coupled_last_filter_test as m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
ROOT=m.ROOT;digest=lambda p:sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=m.STEM);args=ap.parse_args()
source=ROOT/'evidence'/(args.receipt+'.json');raw=json.loads(source.read_bytes());ctx.prec=768
assert raw['status'] in ['complete_last_filter_analytic_estimator_excessive_numerical_stage_not_executed',
 'complete_implemented_last_filter_passes_signed_output_open','complete_implemented_last_filter_estimator_excessive']
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_coupled_last_filter_test.py').read_text(encoding='utf-8')
protected=json.loads(m.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
assert digest(Path(sys.executable))==rt['python_sha256']
assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
r=raw['result'];v={k:dec(x) for k,x in r['constants'].items()};a=r['analytic']
total=lambda obj:sum((dec(x) for x in obj.values()),arb(0))
assert total(a['scalar_residual_debits']).overlaps(v['complete_scalar_residual'])
assert total(a['TT_residual_debits']).overlaps(v['complete_TT_residual'])
assert total(a['quantum_debits']).overlaps(dec(a['TT_residual_debits']['complete_actual_minus_reference']))
assert len(a['quantum_debits'])==10 and len(a['scalar_residual_debits'])==3
qt=arb('.338')/(1-m.cs('current-filtered-mixed-green-residual')['epsilon']);ft=qt/(1+qt)
assert (v['complete_TT_residual']+qt.sqrt()*v['complete_scalar_residual']).overlaps(v['complete_analytic_filter_error'])
assert (ft*ft*arb('1e33')).overlaps(v['last_filter_screen'])
assert (ft*v['implemented_input_energy_error']).overlaps(v['inherited_input_filter_error'])
assert (v['complete_analytic_filter_error']+v['inherited_input_filter_error']).overlaps(v['analytic_plus_implemented_input_error'])
assert all(dec(x)>0 for x in a['filter_polynomial_positivity_guards'])
p,c,bs=m.fr.load();m.lr.model_rows(p,c);c=m.setup(p,c)
assert (v['last_filter_to_original_bolt_multiplier']*ft*ft).overlaps(c['first_residual_to_bolt_multiplier'])

# The ordinary cancellation is polynomial on the entire spectrum. These
# fixtures independently compare the generated rational profiles with it.
def rational(name,x,stage):
 out=arb(0)
 for coeff,poles in (m.old.terms(p,c,name) if stage==2 else m.terms(p,c,name)):
  den=(x*x+c['d1']*x+c['d0'])**(stage-1)
  for mu in poles:den*=x+mu
  out+=m.pv(coeff,x)/den
 return out
fixtures=[]
for x in [arb(0),arb('.1'),arb(1),arb(10),arb('1e6')]:
 F=c['kappa']*(c['Cv']+c['delta']*x)/(x*x+c['d1']*x+c['d0'])
 f2=rational('f',x,2);f3=rational('f',x,3);diff=rational('diff',x,3)
 xi=rational('xidiff',x,3);sig=rational('sig',x,3)
 assert f3.overlaps(F*f2) and (f3+diff).overlaps(f2)
 assert xi.overlaps(-2*(x+2*p['a'])*diff/(3*(x+c['muY'])))
 assert ((x+c['muS'])*sig).overlaps(c['scalar_c']*xi)
 xi3=-2*(x+2*p['a'])*f3/(3*(x+c['muY']))
 tau0=m.cs('current-scalar-green-source-residual')['tau0']
 assert (c['classical_A']*(x+c['muT'])*xi3).overlaps(4*tau0/(p['beta']*p['L']**4)*(c['Cv']+c['delta']*x)*sig)
 fixtures.append(dict(p=enc(x),both_equations=True,full_rational_profile=True))

cells=a['reconstruction_cells'];sc=arb('.0004');touch=0;ex=ew=arb(0)
assert dec(cells[0]['s0']).contains(0) and dec(cells[-1]['s1']).overlaps(p['L'])
for left,right in zip(cells,cells[1:]):assert left['s1']==right['s0']
for cell in cells:
 s0,s1=dec(cell['s0']),dec(cell['s1'])
 if s0.overlaps(sc) or s1.overlaps(sc):touch+=1;assert cell['outside']
 if not(s1<sc):assert cell['outside']
 ex+=dec(cell['xi_diff_error'])**2;ew+=dec(cell['xi_weighted_error'])**2
assert touch==2 and ex.sqrt().overlaps(v['exi']) and ew.sqrt().overlaps(v['exiw'])
assert len(a['spectral_cells'])==512
assert len(a['local_return']['whole_cells'])==152 and len(a['local_return']['outer_rows'])==4
local=a['local_return']['constants']
assert (dec(local['complete_classical_two_row_difference'])+dec(local['complete_local_quantum_rows'])+
 dec(local['complete_adjoint_York_return'])).overlaps(dec(a['TT_residual_debits']['complete_local_return']))

numerical=r['scope']['numerical_last_filter'];product_checks=[]
if numerical:
 assert v['original_X_to_Et_comparison'].overlaps(m.cs('current-filter-numerical-representation-attempt3')['original_X_to_Et_comparison'])
 assert (v['parameter_Et_error']+v['original_X_to_Et_comparison']*(v['contour_X_error']+v['node_and_weight_X_error'])).overlaps(v['numerical_output_error'])
 assert (v['analytic_plus_implemented_input_error']+v['numerical_output_error']).overlaps(v['complete_implemented_last_filter_error'])
 assert (v['last_filter_screen']-v['complete_implemented_last_filter_error']).overlaps(v['remaining_last_filter_account'])
 prior=m.cs('current-filter-numerical-representation-attempt3')
 extra=v['last_filter_to_original_bolt_multiplier']*(v['complete_analytic_filter_error']+v['numerical_output_error'])
 assert extra.overlaps(v['incremental_global_bolt_error_charged'])
 assert (prior['cumulative_known_error_subtotal_NOT_complete_output']+extra).overlaps(v['known_cumulative_error_subtotal_NOT_complete_bolt'])
 assert r['scope']['complete_last_filter_screen_passes']==bool(v['complete_implemented_last_filter_error']<v['last_filter_screen'])
 assert len(r['parameter_error']['whole_spectral_cells'])==1024
 assert len(r['parameter_error']['cross_numerator'])==12
 assert len(r['numerical_field']['upper_half_nodes'])==384
 assert len(r['contour_bounds'])==3
 actual,nominal,ad,nd=m.profile_parameters(p,c)
 assert len(nd[0])==4 and len(nd[1])==9
 Q=lambda q:Fraction(int(q['numerator']),int(q['denominator']))
 field=r['numerical_field'];assert Q(field['constant'])==0
 base=[Q(q) for q in field['base_poles']]
 # Exact rational B4 moments are structural, not tested only at Q samples.
 residues=[1/(m.polyeval(m.polypoles([v for k,v in enumerate(base) if k!=j]),-mu)) for j,mu in enumerate(base)]
 for power in range(3):assert sum(w*mu**power for w,mu in zip(residues,base))==0
 for node in [field['upper_half_nodes'][j] for j in [0,127,128,255,256,383]]:
  mu=m.num.cball((Q(node['pole']['real']),Q(node['pole']['imag'])))
  assert mu.real>0 and mu.imag>0
  poles=[m.acb(m.qball(x)) for x in base]+[mu,m.acb(c['muY'])]
  weights=[]
  for j,z in enumerate(poles):
   den=m.acb(1)
   for k,w in enumerate(poles):
    if k!=j:den*=w-z
   weights.append(1/den)
  for power in range(5):assert sum((w*z**power for w,z in zip(weights,poles)),m.acb(0)).contains(0)
  product_checks.append(dict(positive_real_pole=True,product_moments=True))
 assert len(r['point_and_endpoint_queries'])==5
 for query in r['point_and_endpoint_queries']:
  assert all(dec(x).is_finite() for x in query['values'].values())
 assert all(dec(x).is_zero() for x in r['point_and_endpoint_queries'][-1]['values'].values())
else:
 assert r['numerical_output_error'] is None and r['numerical_field'] is None
 assert not v['analytic_plus_implemented_input_error']<v['last_filter_screen']

assert r['full_signed_Gt'] is None
for key in ['full_signed_response','full_bolt_certificate','original_P32_defect','actual_mass_column',
 'original_graph_physical_boundary_admission','nonlinear_existence','compatibility','physical_matching',
 'validity_noise','all_sector_stability','construction_obstruction','goal_complete']:
 assert r['scope'][key] is False
out=ROOT/'evidence'/(args.receipt+'-verification.json');assert not out.exists()
record=dict(status='complete_last_filter_argument_and_composed_error_accounting_reviewed',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(source),review_source_sha256=digest(Path(__file__)),
 frozen_inputs_checked=len(raw['input_sha256_after']),protected_paths_checked=len(protected),
 review_type='author ordinary argument, both equations, product identities and composed output accounting',
 independent_specialist_review=False,lean_verification=False,full_equation_fixtures=fixtures,
 whole_reconstruction_cells=len(cells),whole_local_cells=152,whole_spectral_cells=512,
 cold_split_touching_cells_charged=touch,numerical_last_filter_reviewed=numerical,product_fixtures=product_checks,
 complete_last_filter_screen_passes=r['scope']['complete_last_filter_screen_passes'],full_signed_response=False,goal_complete=False)
out.write_bytes((json.dumps(record,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in record.items() if not isinstance(v,list)},indent=2))
