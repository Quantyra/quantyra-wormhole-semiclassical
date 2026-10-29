"""Author continuum-accounting, independent reconstruction and provenance review."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
from math import comb
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_source_adapted_signed_response_enclosure as m
arb,ctx,dec,enc,J=m.arb,m.ctx,m.dec,m.enc,m.J
ROOT=m.ROOT;STEM=m.STEM
def digest(p):return sha256(p.read_bytes()).hexdigest()
ctx.prec=512
source=ROOT/'evidence'/(STEM+'-attempt3.json');raw=json.loads(source.read_bytes())
assert raw['status']=='complete_exact_analytic_actual_first_residual_passes_caller_and_later_filters_open'
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_source_adapted_signed_response_enclosure.py').read_text(encoding='utf-8')
protected=json.loads(m.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
assert digest(Path(sys.executable))==rt['python_sha256']
assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
r=raw['result'];v={k:dec(x) for k,x in r['constants'].items()};d={k:dec(x) for k,x in r['debits'].items()}
assert len(d)==10 and 'full_local_annular_kernel' in d and 'actual_constant_contact_OUTSIDE_cold' in d
assert sum(d.values(),arb(0)).overlaps(v['complete_actual_minus_reference_response_upper'])
assert v['complete_actual_first_residual_upper']<arb('8.669e32')
assert v['complete_actual_minus_reference_response_upper']<arb('3.024e32')
assert v['remaining_UNEXECUTED_caller_Et_account']>arb('1.309e32')
assert v['actual_source_L2_upper']<arb('2.838e48') and v['actual_source_gradient_upper']<arb('1.639e48')
assert v['reconstruction_L2_error']<arb('1e22') and v['reconstruction_gradient_error']<arb('1e24')
assert v['complete_outside_cold_source_L2_upper']<arb('1e23')
p,c,bs=m.fr.load();L,a,ga=[p[k] for k in ['L','a','gamma']];sc=arb('.0004')
assert 1/(L-sc)<arb('1.001')
assert (2*(2*ga*sc+a*sc*sc)).sqrt()/2<arb('.01')
assert (arb(2).sqrt()/128).sinh()/(arb(2).sqrt()/128)<2
assert r['cell_count']==len(r['cells']) and dec(r['cells'][0]['s0']).contains(0)
assert dec(r['cells'][-1]['s1']).overlaps(L)
for left,right in zip(r['cells'],r['cells'][1:]):assert left['s1']==right['s0']
touch=0;tail2=arb(0)
for cell in r['cells']:
 s0,s1=dec(cell['s0']),dec(cell['s1'])
 if s1.overlaps(sc) or s0.overlaps(sc):
  touch+=1;assert cell['outside_cold']
 if not (s1<sc):assert cell['outside_cold']
 if cell['outside_cold']:tail2+=dec(cell['reconstruction_error_L2'])**2
assert touch==2
ptail=m.lr.profile_bounds(p,c,m.lr.profile_components(p,c,'f'),sc,L)
assert (ptail+tail2.sqrt()).overlaps(v['complete_outside_cold_source_L2_upper'])
# Independent beta-integral derivation after p=m*t/(1-t).
assert sum((Fraction((-1)**j*comb(4,j),j+1) for j in range(5)),Fraction(0))==Fraction(1,5)
assert sum((Fraction((-1)**j*comb(3,j),j+2) for j in range(4)),Fraction(0))==Fraction(1,20)
# Check the new coefficient assembly against the UNFACTORED accepted field,
# differentiating local ODE jets rather than the new coefficient formula.
fixtures=[]
for xp in [arb('.3'),arb(3),sc/(2*ga),arb('.9')*L/ga]:
 x=J([xp,arb(1)]+[arb(0)]*5);z=L-ga*x;q2=2*x+a*x*x
 q=2*x+a*x*x+sum((bs[j]*ga**(j-2)*x**j for j in range(3,7)),0)
 R=q/q2;Dx=z*q.d().d().d()+3*ga*q.d().d();mu=4*a/3
 A=-(ga*x/L)*sum(((z/L)**j for j in range(5)),0)+(z/L)**5*(R-1)
 B=z**4/L**5*(2*z*(2*q2.d()*R.d()+q2*R.d().d())/3+8*ga*q.d()/3)
 C=z**4/L**5*(5*ga*R-z*R.d());E=-2*z**4/(3*L**5)
 for i in range(5):
  rr,hr,ff,hf,ww=[arb(int(k==i)) for k in range(5)]
  rj=J([rr,hr/q2.c[0]]);hj=J([hr,mu*rr-ff]);fj=J([ff,hf/q2.c[0]])
  wj=J([ww,Dx.c[0]*rr]);MM=z*q.d().d()+4*ga*q.d()
  direct=z**4/L**5*(z*R*fj+(2*MM/3-z*R*mu)*rj+(5*ga*R-z*R.d())*hj-2*wj/3)-fj
  value=B.c[0]*rr+C.c[0]*hr+A.c[0]*ff+E.c[0]*ww
  root=q2.c[0].sqrt()
  gradient=root*(B.d()+mu*C+E*Dx).c[0]*rr+(B+q2*C.d()).c[0]*hr/root
  gradient+=root*(A.d()-C).c[0]*ff+A.c[0]*hf/root+root*E.d().c[0]*ww
  assert value.overlaps(direct.c[0]) and gradient.overlaps(root*direct.d().c[0])
  fixtures.append(dict(x=enc(xp),basis=i,field_agrees=True,gradient_agrees=True))
for key in ['actual_signed_nonzero_center','numerical_caller_error','later_filter_application_errors','signed_Gt']:assert r[key] is None
assert r['proposed_t_5_over_4_trace_theorem_proved'] is False
assert not any(r['scope'][k] for k in ['actual_H5_source_assumed','numerical_first_inverse_applied',
 'full_P32_or_Gt_certified','original_admission','nonlinear_existence','compatibility','physical_matching',
 'validity_noise','all_sector_stability','construction_obstruction','goal_complete'])
dev=[]
for suffix in ['', '-attempt2']:
 path=ROOT/'evidence'/(STEM+suffix+'.json');dr=json.loads(path.read_bytes())
 dev.append(dict(path=str(path.relative_to(ROOT)),sha256=digest(path),status=dr['status'],
  acceptance=False,reason='Superseded by explicit annulus debit and conservative split-touching cell coverage.'))
out=ROOT/'evidence'/(STEM+'-verification.json');assert not out.exists()
review=dict(status='complete_actual_first_analytic_residual_reviewed_later_applications_open',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),review_type='author ordinary-argument and independent reconstruction accounting',
 independent_specialist_review=False,lean_verification=False,receipt_sha256=digest(source),
 review_source_sha256=digest(Path(__file__)),frozen_inputs_checked=len(raw['input_sha256_after']),
 protected_paths_checked=len(protected),whole_cells_checked=len(r['cells']),split_touching_cells_charged=touch,
 independent_reconstruction_fixtures=fixtures,preserved_development=dev,
 continuum_scope=['The existing paired common-regulator form theorem defines the actual fixed-trial response.',
  'The cold cap bound is reused only on its verified region; the whole complementary input tail is charged.',
  'The kernel-only reserve, finite second contact, annulus, both global returns and common local restoration are separate.',
  'Both York equations, all rational profile factors, physical dz and the Volterra return remain.',
  'The actual Et residual equals the energy error for the first actual Riesz inverse; no reference factor is omitted.',
  'No heat-rate theorem, numerical caller error, later filter or original geometric admission is inferred.'])
out.write_bytes((json.dumps(review,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=review['status'],inputs=review['frozen_inputs_checked'],protected=len(protected),
 cells=len(r['cells']),fixtures=len(fixtures),touching_cells=touch),indent=2))
