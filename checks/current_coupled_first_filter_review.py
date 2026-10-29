"""Author argument, independent two-row reconstruction and ledger review."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_coupled_first_filter_test as m
arb,ctx,dec,enc,J=m.arb,m.ctx,m.dec,m.enc,m.J
ROOT=m.ROOT;digest=lambda p:sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=m.STEM);args=ap.parse_args()
source=ROOT/'evidence'/(args.receipt+'.json');raw=json.loads(source.read_bytes());ctx.prec=512
assert raw['status'] in ['complete_analytic_coupled_filter_passes_caller_and_signed_output_open',
                         'complete_analytic_coupled_filter_estimator_excessive']
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_coupled_first_filter_test.py').read_text(encoding='utf-8')
protected=json.loads(m.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
assert digest(Path(sys.executable))==rt['python_sha256']
assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
r=raw['result'];v={k:dec(x) for k,x in r['constants'].items()}
total=lambda name:sum((dec(x) for x in r[name].values()),arb(0))
assert total('scalar_residual_debits').overlaps(v['complete_scalar_residual'])
assert total('TT_residual_debits').overlaps(v['complete_TT_residual'])
assert total('quantum_debits').overlaps(dec(r['TT_residual_debits']['complete_actual_minus_reference']))
assert len(r['quantum_debits'])==10 and len(r['scalar_residual_debits'])==3
qt=arb('.338')/(1-m.cs('current-filtered-mixed-green-residual')['epsilon']);ft=qt/(1+qt)
assert (v['complete_TT_residual']+qt.sqrt()*v['complete_scalar_residual']).overlaps(v['complete_filter_energy_error'])
assert (ft*arb('1e33')).overlaps(v['analytic_filter_screen'])
assert (v['filter_error_to_original_bolt_multiplier']*v['complete_filter_energy_error']).overlaps(v['filter_error_original_bolt_component'])
assert (v['analytic_filter_screen']-v['complete_filter_energy_error']).overlaps(v['remaining_UNEXECUTED_caller_Et_account'])
assert r['scope']['analytic_filter_screen_passes']==bool(v['complete_filter_energy_error']<v['analytic_filter_screen'])
assert all(dec(x)>0 for x in r['filter_polynomial_positivity_guards'])
p,c,bs=m.fr.load();m.lr.model_rows(p,c);c=m.setup(p,c);L,a,ga=[p[k] for k in ['L','a','gamma']]
assert abs(c['h']-arb(4)/3)<arb('.01') and abs(c['e'])<arb('.001')
assert (c['classical_A']*c['kappa']-4*m.cs('current-scalar-green-source-residual')['tau0']*c['scalar_c']/(p['beta']*L**4)).contains(0)

# Independent rational identities at fixtures, supplementing the polynomial
# cancellation and entire-spectrum positivity proof in the ordinary argument.
def rational(name,x):
 out=arb(0)
 for num,poles in m.terms(p,c,name):
  den=x*x+c['d1']*x+c['d0']
  for pole in poles:den*=x+pole
  out+=m.pv(num,x)/den
 return out
fixtures=[]
ll=a*p['ell']*(p['ell']+1);ln=c['large_mass_lambda']
for x in [arb(0),arb('.1'),arb(1),arb(10),arb('1e6')]:
 hfac=1+(c['Cv']-1)*ln/(x+ln)
 f1=c['K0']*hfac/((x+c['muT'])*(x+c['muS'])*(x+ll))
 denom=(x+c['muT'])*(x+c['muS'])+c['kappa']*(c['Cv']+c['delta']*x)
 filt=c['kappa']*(c['Cv']+c['delta']*x)/denom
 f2=rational('f',x);diff=rational('diff',x);xi=rational('xidiff',x);sig=rational('sig',x)
 assert f2.overlaps(filt*f1) and diff.overlaps((1-filt)*f1)
 assert xi.overlaps(-2*(x+2*a)*diff/(3*(x+c['muY'])))
 assert ((x+c['muS'])*sig).overlaps(c['scalar_c']*xi)
 xi2=-2*(x+2*a)*f2/(3*(x+c['muY']))
 tau0=m.cs('current-scalar-green-source-residual')['tau0']
 assert (c['classical_A']*(x+c['muT'])*xi2).overlaps(4*tau0/(p['beta']*L**4)*(c['Cv']+c['delta']*x)*sig)
 fixtures.append(dict(p=enc(x),both_reference_equations=True,complete_high_mass_profile=True))

# Direct canonical nu and QJ/z reconstruction, independently of the new
# factored xi and xi/z^2 comparison coefficients.
xf=[]
for xp in [arb('.3'),arb(3),arb('.0002')/ga,arb('.9')*L/ga]:
 x=J([xp,arb(1)]+[arb(0)]*5);z,q2,R,q=m.fr.geometry_jets(p,bs,x)
 MM=z*q.d().d()+4*ga*q.d()
 A=(z/L)**5*R-1
 B=z**4/L**5*(2*MM/3-z*R*c['muY'])
 C=z**4/L**5*(5*ga*R-z*R.d());E=-2*z**4/(3*L**5)
 for i in range(4):
  rr,hr,ff,ww=[arb(int(k==i)) for k in range(4)]
  nu=z**4/(3*L**5)*(-z*R*ff+(z*R*c['muY']-MM)*rr-(7*ga*R-z*R.d())*hr+ww)
  qJ_over_z=-2*ga*z**4*R*hr/(3*L**5)
  direct=2*nu-2*qJ_over_z;reference=-2*ff/3-4*a*rr/9
  factored=-2*A*ff/3-(B+4*a*A/9)*rr-2*C*hr/3-E*ww
  assert (direct.c[0]-reference).overlaps(factored.c[0])
  weighted=direct/z**2-reference/L**2
  from_difference=factored/z**2+reference*(1/z**2-1/L**2)
  assert weighted.c[0].overlaps(from_difference.c[0])
  xf.append(dict(x=enc(xp),basis=i,xi_and_weighted_xi_agree=True))

cells=r['reconstruction_cells'];sc=arb('.0004')
assert dec(cells[0]['s0']).contains(0) and dec(cells[-1]['s1']).overlaps(L)
for left,right in zip(cells,cells[1:]):assert left['s1']==right['s0']
touch=0;ex=ew=arb(0)
for cell in cells:
 s0,s1=dec(cell['s0']),dec(cell['s1'])
 if s0.overlaps(sc) or s1.overlaps(sc):touch+=1;assert cell['outside']
 if not(s1<sc):assert cell['outside']
 ex+=dec(cell['xi_diff_error'])**2;ew+=dec(cell['xi_weighted_error'])**2
assert touch==2 and ex.sqrt().overlaps(v['exi']) and ew.sqrt().overlaps(v['exiw'])
assert r['local_return']['cell_count']==len(r['local_return']['cells'])
assert len(r['spectral_cells'])==512
for key in ['numerical_caller_error','remaining_filter_error','complete_signed_Gt']:assert r[key] is None
for key in ['numerical_filter_application','full_bolt_certificate','original_P32_defect','actual_mass_column',
 'original_graph_physical_admission','nonlinear_existence','compatibility','physical_matching',
 'validity_noise','all_sector_stability','construction_obstruction','goal_complete']:
 assert r['scope'][key] is False

out=ROOT/'evidence'/(args.receipt+'-verification.json');assert not out.exists()
review=dict(status='complete_coupled_analytic_filter_argument_and_ledger_reviewed',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(source),
 review_source_sha256=digest(Path(__file__)),frozen_inputs_checked=len(raw['input_sha256_after']),
 protected_paths_checked=len(protected),review_type='author ordinary argument and independent reconstruction accounting',
 independent_specialist_review=False,lean_verification=False,
 full_reference_equation_fixtures=fixtures,independent_actual_two_row_fixtures=xf,
 whole_reconstruction_cells=len(cells),whole_local_row_cells=r['local_return']['cell_count'],
 spectral_integral_whole_cells=len(r['spectral_cells']),cold_split_touching_cells_charged=touch,
 analytic_screen_passes=r['scope']['analytic_filter_screen_passes'],numerical_caller_certified=False,
 full_signed_response=False,goal_complete=False)
out.write_bytes((json.dumps(review,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in review.items() if not isinstance(v,list)},indent=2))
