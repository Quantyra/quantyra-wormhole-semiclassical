"""Actual-domain proper-gauge comparison test, not a gravitational solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial, comb
import json, sys, traceback
if '--arb-path' in sys.argv: sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-proper-gauge-neighborhood-test'
OWN=['notes/'+STEM+'.md','checks/current_proper_gauge_neighborhood_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/'+n+'.md' for n in ['stationary-centered-canonical-response',
 'finite-outer-admission-and-response','current-rough-paired-variation',
 'current-original-admission-map-budget-test','current-native-output-structured-source-test',
 'current-schur-residual-phase-domain','current-global-mixed-schur-preconditioner',
 'current-controlled-solution-audit-closeout']]+[
 'evidence/current-global-mass-conversion.json','evidence/current-localized-bolt-adjoint-sensitivity.json']
CHECKS=[]
def sha(p): return sha256(p.read_bytes()).hexdigest()
def check(n,ok):
 CHECKS.append({'name':n,'passed':bool(ok)});print(n+(': passed' if ok else ': FAILED'),flush=True)
 if not ok: raise AssertionError(n)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 if isinstance(x,dict):
  d=x['decimal_enclosure'];return arb(d['mid']+' +/- '+d['rad'])*arb(10)**int(d['exponent'])
 return arb(x.strip('[]'))
def rat(x):return arb(int(x.p))/int(x.q)
def mul(a,b,n):return [sum((a[k]*b[j-k] for k in range(j+1) if k<len(a) and j-k<len(b)),arb(0)) for j in range(n)]
def expjet(a,n):
 out=[a[0].exp()]
 for j in range(1,n):out.append(sum((k*a[k]*out[j-k] for k in range(1,j+1)),arb(0))/j)
 return out
def derivative(a):return [(j+1)*a[j+1] for j in range(len(a)-1)]
def Lpower(w,fac,n):
 out=w[:]
 for _ in range(n):out=mul(fac,derivative(out),len(out)-1)
 return out[0]
def bernstein_guard(poly,t):
 p=sp.Poly(poly,t);n=p.degree()
 return max(abs(sum(p.nth(j)*sp.Rational(comb(k,j),comb(n,j)) for j in range(k+1))) for k in range(n+1))
def execute():
 ctx.prec=512
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_frozen',all(sha(ROOT/p)==h for p,h in rev['scientific_sha256'].items()))
 rr=[json.loads((ROOT/p).read_bytes())['runs'][-1] for p in INPUTS[-2:]]
 check('actual_parameter_receipts_accepted',all(r['status']=='passed' for r in rr))
 pars=rr[0]['result']['parameters'];bp=rr[1]['result']['parameters']
 be,u,ee,tau=[dec(pars[k]) for k in ['beta','u','eps','tau']];h=dec(bp['h']);s0=h/(1+h)
 check('actual_beta_guard',be>arb('.16'))
 qco=mul([arb(1),arb(-2),arb(1)],[be,2*be,3*be+u,arb(0),ee],7);qco[3]-=tau
 derivguards=[]
 for j in range(7):
  derivguards.append(sum((abs(qco[k])*factorial(k)/factorial(k-j)*arb('.87')**(k-j) for k in range(j,7)),arb(0)))
 check('whole_actual_polynomial_derivative_guards',all(x<arb('1e6') for x in derivguards))
 check('proper_cap_stays_in_seed_chart',1-(-arb('2.01')).exp()<arb('.87'))
 check('actual_cap_Q_lower',be*(arb('.13')**2-s0**2)>arb('.002'))
 check('actual_Q_upper',be+u+abs(ee)<arb('.17'))
 qstarlo=be*((1-arb(1)/40)**2-s0**2)
 check('actual_observation_Q_lower',qstarlo>arb('.15'))
 check('actual_observation_proper_interval',arb(1)/40/arb('.17').sqrt()>arb('.06') and arb(1)/40/arb('.15').sqrt()<arb('.063'))
 rz=arb('1e-10');ry=arb('1e-12');qvar=arb('1e6')*((rz).exp()-1)
 check('complex_seed_Q_variation',qvar<arb('.000101'))
 check('seed_positive_branch_preserved',qvar<arb('.002')/2)
 check('seed_proper_flow_contraction',ry*2/rz<arb('.03'))
 check('actual_seed_seven_jet_guard',max(10*factorial(j)*arb('1e12')**j for j in range(1,8))<arb('1e100'))
 # Whole-interval cutoff and direction jets, with exact rational Bernstein coefficients.
 t=sp.Symbol('t');p=4**12*t**12*(1-t)**12
 beta_norm=sp.factorial(23)/(sp.factorial(11)**2)
 step=sp.integrate(beta_norm*t**11*(1-t)**11,t)
 check('cutoff_exact_endpoint_values',step.subs(t,0)==0 and step.subs(t,1)==1)
 check('original_endpoint_jet_matching',all(sp.diff(p,t,j).subs(t,x)==0 for j in range(12) for x in [0,1]))
 cb=[arb(1)]+[rat(bernstein_guard(sp.diff(step,t,j),t))/arb('.005')**j for j in range(1,8)]
 pb=[rat(bernstein_guard(sp.diff(p,t,j),t)) for j in range(6)]
 N=arb('1e200');amp=arb('1e-20');av=arb('1e-12');width=arb('.01')
 Tguards=[amp*N**(j-5)*sum((comb(j,k)*cb[k]/N**k for k in range(j+1)),arb(0)) for j in range(8)]
 check('whole_packet_jet_guard',all(Tguards[j]<2*amp*N**(j-5) for j in range(8)))
 def denom_lower(j,lo,hi):return factorial(j)*arb(100)**j*(lo**(4-j) if j<=4 else 1/hi)
 U=max(av*pb[j]/width**j/denom_lower(j,arb('.01'),arb('.095')) for j in range(6))
 U=max(U,(av*pb[0]+av*pb[1]/width)/100)
 Tin=max(Tguards[j]/denom_lower(j,arb('.045'),arb('.080')) for j in range(6))
 check('canonical_packet_five_jet_admitted',Tin<arb('1e-14'))
 J=av*width*arb(4)**12*factorial(12)**2/factorial(25)
 J2=av**2*width*arb(16)**12*factorial(24)**2/factorial(49)
 check('exact_shift_moments_positive',J>0 and J2>0)
 check('first_moment_integral',sp.integrate(p,(t,0,1))==sp.Rational(4**12*factorial(12)**2,factorial(25)))
 check('second_moment_integral',sp.integrate(p*p,(t,0,1))==sp.Rational(16**12*factorial(24)**2,factorial(49)))
 def jets(kind):
  trig=[1,0,-1,0,1,0,-1,0] if kind=='cos' else [0,1,0,-1,0,1,0,-1]
  T=[amp*N**(j-5)*trig[j]/factorial(j) for j in range(8)]
  factor=expjet([-x for x in T],8)
  W=[-T[0]]+[arb('0 +/- 1e100')/factorial(j)-T[j] for j in range(1,8)]
  return Lpower(W,factor,6),Lpower(W,factor,7)
 cos6,_=jets('cos');sin6,sin7=jets('sin')
 first=abs(-J*cos6);second=abs(J*J*sin7-J2*sin6)
 den=arb('2e14')
 check('original_fifth_weight_upper',120*arb(100)**5*(1-(arb('.04')).log())/arb('.04')<den)
 firstout=first/den;secondout=second/den
 check('finite_first_output_failure',firstout>arb('1e100'))
 check('finite_second_output_failure',secondout>arb('1e250'))
 check('finite_normalized_first_failure',firstout/U>arb('1e100'))
 check('finite_normalized_second_failure',secondout/U**2>arb('1e250'))
 nonlinear_error=arb('1e120')*amp/N
 check('whole_finite_endpoint_chain_remainder',nonlinear_error<amp)
 background=3*amp/arb('.04')**4
 check('actual_proper_background_admitted',background<arb('2e-14') and background<arb('1.5e-7'))
 # Exact moving-coordinate first and second derivatives by implicit differentiation.
 hh,ju,jv,juv,tp,wp,wpp,uu,vv,up,vp=sp.symbols('h Ju Jv Juv Tp Wp Wpp u v up vp')
 yu=-hh*ju;yv=-hh*jv
 yuv=-tp*yu*yv-uu*yv-vv*yu-hh*juv
 full=-up*yv-vp*yu+wpp*yu*yv+wp*yuv
 expected=hh*(up*jv+vp*ju)+hh**2*(wpp-tp*wp)*ju*jv+hh*wp*(uu*jv+vv*ju-juv)
 check('complete_coordinate_second_variation',sp.expand(full-expected)==0)
 check('separated_direction_second_contact',sp.expand(full.subs({uu:0,vv:0,up:0,vp:0,ju:sp.Symbol('J'),jv:sp.Symbol('J'),juv:sp.Symbol('J2')})-(hh**2*(wpp-tp*wp)*sp.Symbol('J')**2-hh*wp*sp.Symbol('J2')))==0)
 xx=sp.Symbol('x',positive=True)
 check('log_embedding_integral',sp.integrate(1/(xx*(1+sp.log(xx))**2),(xx,1,sp.oo))==1)
 return {'precision_bits':512,'parameters':{'N':enc(N),'amplitude':enc(amp)},
  'finite_actual_domain_test':{k:enc(v) for k,v in {'canonical_packet_U5_upper':Tin,
   'direction_U5_upper':U,'J':J,'J2':J2,'proper_background_Nout_upper':background,
   'first_output_weighted_lower':firstout,'second_output_weighted_lower':secondout,
   'normalized_first_derivative_lower':firstout/U,'normalized_second_derivative_lower':secondout/U**2,
   'whole_endpoint_chain_remainder':nonlinear_error,'H8log_to_C7_constant':1/arb.pi().sqrt()}.items()},
  'scope':{'uniform_moving_gauge_C5_derivative_bound_exists':False,
   'original_fixed_endpoint_theorem_rejected':False,'geometric_five_jet_admission_rejected':False,
   'structured_F3_map_rejected':False,'polynomial_higher_jet_constants':None,
   'full_fixed_coordinate_quantum_neighborhood_bound':None,'complete_coupled_bound_certified':False,
   'actual_quantum_source_or_polynomial_application_executed':False,'new_physical_gate_closed':False,
   'construction_obstruction':False,'full_goal_complete':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sha(ROOT/p) for p in INPUTS}
 run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
  'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sha(Path(sys.executable)),
             'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
  check('inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
 out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
 print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
