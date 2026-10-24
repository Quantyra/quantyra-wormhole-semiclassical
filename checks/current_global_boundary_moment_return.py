"""Original-row Mellin extraction and conditional strong return targets.

The actual global remainder is NOT computed by this checker.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
from flint import arb,ctx
import sympy as sp
import current_product_matched_source_projection as old
STEM='current-global-boundary-moment-return'
OWN=['notes/'+STEM+'.md','checks/current_global_boundary_moment_return.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded

def check(name,ok,detail=''):
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rat(x):return arb(int(sp.numer(x)))/int(sp.denom(x))

def algebra():
 u,s,e,g0,b0,bb,der=sp.symbols('u s e g0 B0 B4 Bprime')
 S=35*u**4-84*u**5+70*u**6-20*u**7;d=sp.Rational(1,200)
 check('source_split_endpoints',S.subs(u,0)==0 and S.subs(u,1)==1)
 check('source_split_positive_derivative',sp.simplify(sp.diff(S,u)-140*u**3*(1-u)**3)==0)
 check('source_split_C3',all(sp.diff(S,u,j).subs(u,x)==0 for j in [1,2,3] for x in [0,1]))
 integrals=[];P=S
 for j in range(4):
  if j:P=sp.expand(sp.diff(P,u)-d*P)
  value=sp.integrate(P*P,(u,0,1));integrals.append(value)
  check('positive_exact_transition_integral_'+str(j),value>0)
 cstar=-(-g0/e+sp.Symbol('M'))/der
 delta=sp.simplify((b0+cstar).subs(g0,bb*b0))
 check('original_moment_residue_formula',sp.simplify(delta-(b0*(1+bb/(e*der))-sp.Symbol('M')/der))==0)
 # A quadratic polynomial independently verifies the integral sign.
 t,A,B,C=sp.symbols('t A B C');poly=A+B*t+C*t*t
 identity=sp.integrate(t*sp.diff(poly,t,2),(t,0,e))
 check('stable_pole_cancellation_identity',sp.simplify(poly.subs(t,0)+e*sp.diff(poly,t).subs(t,e)-poly.subs(t,e)-identity)==0)
 # Symbolic Fourier derivative ledger for the unchanged phase norm.
 p=sp.Integer(1);jet=[]
 for j in range(6):
  if j:p*=s-(j-1)
  translated=sp.Poly(sp.expand(p.subs(s,4+d+u)),u)
  bound=sum(abs(v) for v in translated.all_coeffs())
  check('phase_Sobolev_polynomial_'+str(j),bound<=factorial(j)*100**j)
  jet.append(str(bound))
 # Whole-strip/contour arithmetic from the ordinary complex estimates.
 r=sp.Rational(1,100)
 check('low_strip_bound',21*(100+18)<3000)
 check('shifted_contour_low_bound',21*(200+18)<5000)
 check('pole_cleared_second_derivative',2*(1+r*17)+2*(5+2*r)*(17+r*4)+(6+5*r+r*r)*(8+r*4)<300)
 check('high_frequency_real_part',5*(6+5*r+r*r)+(5+2*r)*5<61)
 check('root_second_to_first_derivative', (26+2*(5+2*r)*(sp.Integer(3)*10**34+120)+(6+5*r+r*r)*300)/(4*(sp.Integer(3)*10**34-120))<4)
 SYMBOLIC.update(transition_integrals=[str(v) for v in integrals],phase_polynomial_sums=jet,
  actual_global_moment_formula='Delta=B0[1+B(4)/(e_* Bprime)]-M_G/Bprime',
  stable_mismatch_bound='2 abs(B0) e_*+abs(M_G)/Bprime',
  original_global_source_terms_retained=['seed source','mass column','actual-minus-normal kernel','full cutoff commutator'])
 return integrals

def execute():
 integrals=algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 prior=json.loads((ROOT/'evidence/current-coupled-outer-return-test.json').read_bytes())['runs'][-1]
 check('prior_full_normal_test_passed',prior['status']=='passed')
 rr=prior['result'];e=old.decode(rr['exact_root_enclosure']);B0=old.decode(rr['full_normal_shape_particular'])
 C,B2,be,la,Fl=p['C'],p['B2'],p['beta'],p['g'],p['Fl'];K=C*be*be/(15*arb.pi()*B2)
 Fs=Fl+23*la/(180*arb.pi());L=10*arb.pi()*Fs/la;d0=1/(480*arb.pi()**2)
 AT=L*d0;D0=-120*arb.pi()**2*K;d=arb(1)/200
 qstar=(4+e)*(1+e);Bprime=4*D0*d0*qstar*(-L+arb(0,120))
 check('global_moment_response_positive',Bprime>arb('.444443') and Bprime<arb('.444445'))
 c3=4*(-D0)*AT*d*(1-5000/L)
 c4low=c3/arb(17).sqrt();c4high=8*(-D0)*d0*(arb(16)/17)**2*(2*L).sqrt()
 check('fourth_gain_high_frequency_dominates',c4high<c4low)
 C3=1/c3;C4=1/c4high
 check('entire_frequency_third_gain',C3<1801)
 check('extra_derivative_changes_scale',C4/C3>arb('1e14'))
 a=arb('.1');La=-a.log();factor=a**d/(arb(2).sqrt()*(1+La))
 # g0=B(4) B0 exactly; parameter intervals retain the same action.
 g0=-24*K*B0
 norms={}
 for m in [2,3]:
  finite=sum((comb(m,j)*rat(integrals[j]) for j in range(m+1)),arb(0))
  norms[m]=g0.abs_upper()*a**(-d)*((2*d).exp()*finite+(1+d*d)**m/(2*d)).sqrt()
 knownmoment=g0.abs_upper()*(1+La)*(e*(1+La)).exp()
 co=1/(2*(d-e)).sqrt();cf=1/(2*(arb('.5')+e)).sqrt();pole=B0.abs_upper()*e.abs_upper()
 E3=arb('1e-11');Efar=arb('1e-10');moment=knownmoment+co*E3+cf*Efar
 mismatch=2*pole+moment/Bprime
 remainder3=factor*C3*(E3+norms[3]);total=pole+mismatch/(1+La)+remainder3
 check('conditional_shape_phase_target',total<arb('1e-8'))
 check('conditional_mismatch_target',mismatch<arb('1e-8'))
 check('known_endpoint_H2_remainder_passes',factor*C4*norms[2]<arb('1e-9'))
 half=arb('5e-9');target2=half/(factor*C4)-norms[2];target3=half/(factor*C3)-norms[3]
 check('both_residual_allocations_positive',target2>0 and target3>0)
 return dict(sigma='4.005',exact_slow_exponent_offset=enc(e),Bprime_at_slow_root=enc(Bprime),moment_response=enc(1/Bprime),
  C3_H3_to_H6=enc(C3),C4_H2_to_H6=enc(C4),extra_derivative_improvement_factor=enc(C4/C3),
  canonical_phase_Sobolev_factor=enc(factor),g0=enc(g0),known_source_H2_upper=enc(norms[2]),known_source_H3_upper=enc(norms[3]),
  known_regular_moment_upper=enc(knownmoment),known_H2_phase_remainder_upper=enc(factor*C4*norms[2]),
  outer_moment_weight_factor=enc(co),far_moment_weight_factor=enc(cf),normal_pole_mismatch_upper=enc(2*pole),
  conditional_E3='<=1e-11',conditional_Efar='<=1e-10',conditional_mismatch_upper=enc(mismatch),
  conditional_remainder_phase_upper=enc(remainder3),conditional_total_shape_phase_upper=enc(total),
  half_budget_E2_threshold=enc(target2),half_budget_E3_threshold=enc(target3),
  actual_global_coefficient_well_defined_at_fixed_seed=True,actual_global_coefficient_evaluated=False,
  actual_E3_upper_known=False,actual_Efar_upper_known=False,actual_global_strong_remainder_certified=False,
  source_and_kernel_global_terms_replaced=False,finite_physical_modes_substituted=False,
  conditional_phase_target_is_full_admission=False,actual_controlled_solution=False,
  nonlinear_admission=False,exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False)

def main():
 ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   try:p=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if p.endswith('.py'):paths.add(p)
 before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['symbolic']=SYMBOLIC;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
