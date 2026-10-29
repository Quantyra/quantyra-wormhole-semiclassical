"""Author review of the full representation error, domain and return accounting."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_filter_numerical_representation as m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
ROOT=m.ROOT;digest=lambda p:sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=m.STEM+'-attempt3');args=ap.parse_args()
source=ROOT/'evidence'/(args.receipt+'.json');raw=json.loads(source.read_bytes());ctx.prec=768
assert raw['status'] in ['complete_finite_product_filter_caller_passes','complete_finite_product_filter_caller_estimator_excessive']
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
for key,path in [('ordinary_argument_utf8','notes/'+m.STEM+'.md'),
                 ('arithmetic_source_utf8','checks/current_filter_numerical_representation.py'),
                 ('numerical_operator_source_utf8','operators/paired_filter_profile.py')]:
 assert raw[key]==(ROOT/path).read_text(encoding='utf-8')
protected=json.loads(m.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
assert digest(Path(sys.executable))==rt['python_sha256']
assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
r=raw['result'];v={k:dec(x) for k,x in r['constants'].items()}
CE=v['original_X_to_Et_comparison'];X=sum((v[k] for k in ['parameter_X_error','contour_X_error','coefficient_and_node_X_error']),arb(0))
assert (CE*X).overlaps(v['generic_caller_Et_upper'])
assert (v['direct_parameter_Et_error']+CE*(v['contour_X_error']+v['coefficient_and_node_X_error'])).overlaps(v['complete_numerical_caller_Et_error'])
oc=r['output_comparison'];cc={k:dec(x) for k,x in oc['constants'].items()}
assert cc['energy_feedback']<1
assert ((cc['complete_local_dual']+cc['complete_geodesic_dual']+
 sum((dec(x) for x in oc['quantum_entries'].values()),arb(0)))/(1-cc['energy_feedback'])).overlaps(CE)
assert len(oc['quantum_entries'])==7
assert len(oc['whole_cell_coefficients'])==152
assert len(r['parameter_error']['whole_spectral_cells'])==1024
assert len(r['contour_bounds'])==3 and len(r['numerical_field']['upper_half_nodes'])==384
direct=r['direct_parameter_output']
assert len(direct['whole_cells'])==152 and len(direct['outer_rows'])==4
assert ((dec(direct['complete_local_dual'])+dec(direct['complete_geodesic_dual'])+
 dec(direct['complete_quantum_difference_except_feedback']))/(1-cc['energy_feedback'])).overlaps(v['direct_parameter_Et_error'])
assert v['complete_numerical_caller_Et_error']<arb('3.544e27')
assert v['caller_original_bolt_error']<arb('4.490e-16')
assert v['generic_report_over_target']>2000 and v['generic_to_direct_report_gain']>arb('9.85e6')
old=m.cs('current-coupled-first-filter-test-attempt2')
assert (old['complete_filter_energy_error']+v['complete_numerical_caller_Et_error']).overlaps(v['complete_analytic_plus_caller_filter_error'])
assert (old['filter_error_to_original_bolt_multiplier']*v['complete_numerical_caller_Et_error']).overlaps(v['caller_original_bolt_error'])
passed=v['complete_numerical_caller_Et_error']<v['caller_target']
assert r['scope']['numerical_caller_target_passes']==bool(passed)
assert v['caller_target'].contains(arb('1.70e31'))

# Check the independent product partial-fraction identities, including
# the moments which protect the bolt. These are exact polynomial identities;
# ball fixtures supplement the ordinary whole-spectrum proof.
p,c,bs=m.fr.load();m.lr.model_rows(p,c);c=m.old.setup(p,c)
field=r['numerical_field'];model=m.PairedFilterProfile(p,bs,field)
q0,_,_=model.kernel(m.acb(0),arb(3));q1,_,_=model.kernel(2*p['a'],arb(3))
assert q0.overlaps(arb(2).log()/2) and q1.overlaps(3*q0-1)
def Q(q):return Fraction(int(q['numerator']),int(q['denominator']))
base=[m.qball(Q(q)) for q in field['base_poles']]
assert all(b>0 for b in base)
fixtures=[]
for node in [field['upper_half_nodes'][i] for i in [0,127,128,255,256,383]]:
 pole=m.acb(m.qball(Q(node['pole']['real'])),m.qball(Q(node['pole']['imag'])))
 assert pole.real>0 and pole.imag>0
 poles=[m.acb(x) for x in base]+[pole,m.acb(4*p['a']/3)]
 residues=[]
 for j,mu in enumerate(poles):
  den=m.acb(1)
  for k,other in enumerate(poles):
   if j!=k:den*=other-mu
  residues.append(1/den)
 for power in range(5):assert sum((a*mu**power for a,mu in zip(residues,poles)),m.acb(0)).contains(0)
 for z in [m.acb('.2'),m.acb(2),m.acb(50)]:
  product=m.acb(1)
  for mu in poles:product/=z+mu
  partial=sum((a/(z+mu) for a,mu in zip(residues,poles)),m.acb(0))
  assert partial.overlaps(product)
 fixtures.append(dict(pole=node['pole'],first_five_moments_vanish=True,partial_fraction_identity=True))

# The primitive follows from the two Legendre ODEs. Verify the polynomial
# half independently and check its derivative identity using the Q ODE.
from paired_filter_profile import polynomials,polyval
P=polynomials();primitive_checks=[]
for j,Pj in enumerate(P):
 dP=[k*Pj[k] for k in range(1,len(Pj))];ddP=[k*dP[k] for k in range(1,len(dP))]
 ode=m.polyadd(m.polymul([-1,0,1],ddP or [0]),m.polymul([0,2],dP or [0]),[-j*(j+1)*x for x in Pj])
 assert all(arb(x).contains(0) for x in ode)
 mu=m.acb('.7','.1');w=arb('2.3');q,flux,k=model.kernel(mu,w)
 qp=flux/(w*w-1);qpp=(mu/p['a']*q-2*w*qp)/(w*w-1)
 derivative=(2*w*(polyval(Pj,w)*qp-polyval(dP,w)*q)
 +(w*w-1)*(polyval(Pj,w)*qpp-polyval(ddP,w)*q))/(mu/p['a']-j*(j+1))
 assert derivative.overlaps(polyval(Pj,w)*q)
 primitive_checks.append(dict(degree=j,Legendre_polynomial_ODE=True,return_primitive_derivative=True))

assert [q['label'] for q in r['point_and_endpoint_queries']]==['bolt','thermal_scale','scalar_cutoff_scale','cold_split','original_outer_boundary']
for query in r['point_and_endpoint_queries']:
 assert all(dec(x).is_finite() for x in query['values'].values())
assert all(dec(x).is_zero() for x in r['point_and_endpoint_queries'][-1]['values'].values())
assert r['full_signed_Gt'] is None and r['remaining_filter_error'] is None
for key in ['arbitrary_interpolant_certified','full_bolt_certificate','original_P32_defect','mass_column',
 'original_physical_graph_admission','nonlinear_existence','compatibility','physical_matching',
 'validity_noise','all_sector_stability','construction_obstruction','goal_complete']:
 assert r['scope'][key] is False
development=[]
for file in sorted((ROOT/'evidence').glob(m.STEM+'*.json')):
 if file==source or 'verification' in file.name:continue
 prior=json.loads(file.read_bytes())
 development.append(dict(path=str(file.relative_to(ROOT)),sha256=digest(file),status=prior.get('status'),accepted=False))
out=ROOT/'evidence'/(args.receipt+'-verification.json');assert not out.exists()
record=dict(status='complete_finite_product_representation_argument_and_error_accounting_reviewed',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(source),
 review_source_sha256=digest(Path(__file__)),frozen_inputs_checked=len(raw['input_sha256_after']),
 protected_paths_checked=len(protected),review_type='author ordinary argument and independent product/return accounting',
 independent_specialist_review=False,lean_verification=False,whole_local_cells=152,whole_outer_arcs=64,
 whole_parameter_spectral_cells=1024,contour_nodes=768,product_moment_fixtures=fixtures,
 independent_return_primitive_checks=primitive_checks,preserved_development=development,
 caller_target_passes=bool(passed),complete_signed_response=False,goal_complete=False)
out.write_bytes((json.dumps(record,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in record.items() if not isinstance(v,list)},indent=2))
