"""Selected geometric heat coefficient; no complete quantum projection."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import sys,json,traceback,heapq
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
from flint import arb,ctx
import current_paired_regulator_domain_audit as previous
geom=previous.geom;prior=geom.prior;gy=prior.scalar.prev.york
root=Path(__file__).resolve().parents[1]
STEM='current-smeared-heat-projection-domain';PREVIOUS=previous.STEM
OWN=['notes/'+STEM+'.md','checks/current_smeared_heat_projection_domain.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=geom.encoded
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
gy.check=check;prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()

def strong_kernel(Q,z):
 R=-sp.diff(Q,z,2)
 D=lambda f:sp.diff(Q*sp.diff(f,z),z)
 L=lambda f:D(f)+R*f
 v=-z**4
 out=-(L(D(L(v)))+D(R)*L(v)+D(v*D(R)))/105
 out-=(16*R**3*v+24*R**2*D(v)+24*D(R**2*v)+24*D(R*D(v)))/1008
 return sp.expand(out+sp.Rational(4,63)*v)

def algebra():
 t,lam,b=sp.symbols('t lam b',positive=True)
 k=sp.exp(-t*lam*sp.exp(-2*b))
 check('constant_conformal_smeared_heat_identity',sp.simplify(sp.diff(k,b)+2*t*sp.diff(k,t))==0)
 check('constant_second_vertex_retained',sp.simplify(sp.diff(k,b,2)-(-2*t*sp.diff(sp.diff(k,b),t)))==0)
 cgrad=-sp.Rational(1,140)+sp.Rational(1,480)+sp.Rational(1,160)
 ccubic=sp.Rational(1,315)-sp.Rational(1,240)
 check('a6_gradient_contractions',cgrad==sp.Rational(1,840))
 check('a6_cubic_contractions',ccubic==-sp.Rational(1,1008))
 s=sp.Integer(8)
 moment=lambda nu:2**(s-3)*sp.gamma(s/2)**2*sp.gamma(s/2+nu)*sp.gamma(s/2-nu)/sp.gamma(s)
 I6=(moment(0)-moment(1))/(2*sp.pi**2)
 check('flat_sixth_Bessel_moment',sp.simplify(I6+48/(35*sp.pi**2))==0)
 check('flat_a6_gradient_normalization',sp.simplify(-I6/720*(5*sp.pi/8)-1/(840*sp.pi))==0)
 z=sp.symbols('z',positive=True);Q=sp.Function('Q')(z);h=sp.Function('h')(z);v=-z**4
 R=-sp.diff(Q,z,2);D=lambda f:sp.diff(Q*sp.diff(f,z),z)
 a,c=sp.symbols('a c');u=a*h+c*v;Sb=R-2*D(u)
 density=Q*sp.diff(sp.exp(-2*u)*Sb,z)**2/840-sp.exp(-4*u)*Sb**3/1008
 H=sp.diff(density,a,c).subs({a:0,c:0}).expand()
 EL=sum((-1)**j*sp.diff(sp.diff(H,sp.diff(h,z,j)),z,j) for j in range(4))
 expected=strong_kernel(Q,z)-sp.Rational(4,63)*v
 check('full_variation_matches_strong_kernel',sp.expand(EL-expected)==0)
 product=-(-2)**3*sp.exp(-4*b)/1008
 gp=sp.diff(product,b).subs(b,0)
 check('product_a6_gradient',gp==-sp.Rational(2,63))
 check('product_contact_retained',-2*gp==sp.Rational(4,63))
 check('mixed_boundary_jet_weight_exceeds_a6_a7',3+4>5 and 3+4>6)
 for j in range(3):
  J=((1-z**(j+1))*sp.log(1-z)+sum(z**i/sp.Integer(i) for i in range(1,j+2)))/(j+1)-z**(j+2)/sp.Integer(j+2)
  check('log_weight_primitive_'+str(j),sp.simplify(sp.diff(J,z)-z**j*(-sp.log(1-z)-z))==0)
 aa=sp.Function('a')(z);w=sp.Function('w')(z);O=sp.Function('O')(z)
 hb=aa/2-z*w-z*z*sp.diff(w,z)/2
 check('constraint_pairing_keeps_endpoint',sp.simplify(hb*O-aa*O/2-z*z*sp.diff(O,z)*w/2+sp.diff(z*z*O*w/2,z))==0)
 NUMERICS.append(dict(kind='local_coefficient_argument',gradient_coefficient='1/840',cubic_coefficient='-1/1008',
  product_gradient='-2/63',mixed_a6_boundary_jet_weight='7 > 5',mixed_a7_boundary_jet_weight='7 > 6',
  selected_pair_regulator_argument='Ordinary note sections 1-2, not verified by symbolic checks.',
  C2_action_on_rough_base_family_proved=False))

def J(k,z):
 return ((1-z**(k+1))*(1-z).log()+sum((z**i/i for i in range(1,k+2)),arb(0)))/(k+1)-z**(k+2)/(k+2)
def abrange(coeff):
 lo=min(x.lower() for x in coeff);hi=max(x.upper() for x in coeff)
 upper=max(abs(lo),abs(hi))
 lower=arb(0) if lo<=0<=hi else min(abs(lo),abs(hi))
 return lower,upper

def execute():
 p,h,zh,env=gy.actual_parameters();NUMERICS.extend(gy.NUMERICS)
 z=gy.z;O=strong_kernel(gy.Q,z);Op=sp.diff(O,z)
 coeff=lambda f:[gy.evaluate_polynomial(c,env) for c in reversed(sp.Poly(f,z).all_coeffs())]
 co=coeff(O);cop=coeff(Op);zc=1-arb(2)**-49
 A=arb('3.768e-12');al=arb('4e-28');f0=arb('4.34e-11');sr=arb('3.33e-10')*h
 def cell(lo,hi,bb):
  w0=A*(J(0,hi)-J(0,lo))/2+al*(hi-lo)/2
  w1=A*(J(1,hi)-J(2,hi)-J(1,lo)+J(2,lo))
  w1+=(f0+al)*((hi*hi-lo*lo)/2-(hi**3-lo**3)/3)+sr*(hi**3-lo**3)/6
  assert w0>=0 and w1>=0
  ranges=[abrange(x) for x in bb]
  low=(w0*ranges[0][0]+w1*ranges[1][0])/(4*arb.pi())
  high=(w0*ranges[0][1]+w1*ranges[1][1])/(4*arb.pi())
  return dict(lo=lo,hi=hi,bern=bb,weights=[w0,w1],ranges=ranges,lower=low,upper=high,priority=float((high-low).upper()))
 heap=[];seq=0;lower=arb(0);upper=arb(0)
 def add(item):
  nonlocal seq,lower,upper
  seq+=1;heapq.heappush(heap,(-item['priority'],seq,item));lower+=item['lower'];upper+=item['upper']
 for i in range(256):
  lo=zc*i/256;hi=zc*(i+1)/256
  add(cell(lo,hi,[prior.bern_on_cell(c,lo,hi) for c in [co,cop]]))
 while upper.upper()-lower.lower()>arb('1e-12'):
  assert len(heap)<20000,'adaptive_cell_limit'
  _,_,item=heapq.heappop(heap);lower-=item['lower'];upper-=item['upper']
  halves=[prior.scalar.split(a) for a in item['bern']];mid=(item['lo']+item['hi'])/2
  add(cell(item['lo'],mid,[a[0] for a in halves]));add(cell(mid,item['hi'],[a[1] for a in halves]))
  if len(heap)%1000==0:print(json.dumps(dict(event='a6_envelope_refinement',cells=len(heap),width=str(upper-lower))),flush=True)
 cells=sorted((item for _,_,item in heap),key=lambda c:float(c['lo']))
 check('complete_bulk_coverage',cells[0]['lo']==0 and cells[-1]['hi']==zc and all(a['hi']==b['lo'] for a,b in zip(cells,cells[1:])))
 check('envelope_quadrature_width',upper.upper()-lower.lower()<arb('1e-12'))
 oend=gy.evaluate_polynomial(O,{**env,z:zc}).abs_upper()
 wend=sr+2*(1-zc)/zc*(A*(-(1-zc).log()-zc)+f0+al)
 endpoint=zc*zc*oend*wend/(8*arb.pi())
 otail=abrange(prior.bern_on_cell(co,zc,arb(1)))[1]
 tip=arb('1.2e-8')*(1-zc)*otail/(4*arb.pi())
 U=upper+endpoint+tip
 Fl=1-384*arb.pi()*arb('1e22')*p['g'];Dtarget=arb('2.5e-43')*p['B2']*Fl/(4*p['C']);t0=arb('.001')
 capacity=2*(Dtarget-U*t0)/(t0*t0)
 check('bulk_radial_hypothesis',(1-zc)/zc>arb('1e-15'))
 check('local_coefficient_end_terms_retained',endpoint>0 and tip>0)
 check('actual_geometric_coefficients_retained',p['tau']>0 and p['u']>0 and p['eps']<0)
 NUMERICS.append(dict(kind='a6_actual_pair_envelope',O_polynomial=str(O),O_coefficients=[encoded(x) for x in co],
  Oprime_coefficients=[encoded(x) for x in cop],cold_cut=encoded(zc),B_squared=encoded(p['B2']),C=encoded(p['C']),F_lambda=encoded(Fl),
  h_bolt=encoded(h),s0=encoded(sr),bulk_envelope_lower=encoded(lower),bulk_envelope_upper=encoded(upper),
  endpoint_allowance=encoded(endpoint),omitted_interval_allowance=encoded(tip),U6=encoded(U),
  Dtarget=encoded(Dtarget),t0='.001',leading_cap_allowance=encoded(4*p['C']*U*t0/(p['B2']*Fl)),
  required_B_heat_upper=encoded(capacity),capacity_positive=bool(capacity>0),B_heat_known=False,
  full_heat_cap_certified=False,complete_projection_evaluated=False,
  cells=[dict(lo=encoded(c['lo']),hi=encoded(c['hi']),weights=[encoded(x) for x in c['weights']],
   O_abs_range=[encoded(x) for x in c['ranges'][0]],Oprime_abs_range=[encoded(x) for x in c['ranges'][1]],
   lower=encoded(c['lower']),upper=encoded(c['upper'])) for c in cells]))
 print(json.dumps(dict(event='a6_capacity_test',cells=len(cells),U6=str(U),B_heat_capacity=str(capacity),
  capacity_positive=bool(capacity>0),complete_projection_evaluated=False)),flush=True)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,dependencies=deps,invocation=sys.argv,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 try:
  check('previous_scoped_audit_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',prior.scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  continuum_argument_scope='Ordinary selected-pair argument in note; checks verify algebra and coefficient enclosure only.',
  C2_action_on_rough_base_family_proved=False,complete_signed_projection_evaluated=False,full_heat_cap_certified=False,
  actual_source_failure_proved=False,total_source_improvement_certified=False,nonlinear_solution_constructed=False,
  physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
