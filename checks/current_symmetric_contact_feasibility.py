"""Actual probe degeneracy and conditional symmetric-contact graph thresholds."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-symmetric-contact-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_symmetric_contact_feasibility.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def evaluate(poly,x):
 value=arb(0)
 for c in reversed(poly):value=value*x+c
 return value
def derivative(poly):return [i*poly[i] for i in range(1,len(poly))]
def bernstein(poly,a):
 n=len(poly)-1;power=[v*a**j for j,v in enumerate(poly)]
 return [sum((power[j]*comb(k,j)/comb(n,j) for j in range(k+1)),arb(0)) for k in range(n+1)]

def execute():
 pars,h,zh,env=old.scalar.prev.york.actual_parameters();a=arb('.01')
 be,u,ep,ta=pars['beta'],pars['u'],pars['eps'],pars['tau']
 Q=[be,arb(0),u,-4*be-2*u-ta,3*be+u+ep,-2*ep,ep]
 Q1=derivative(Q);Q2=derivative(Q1);Q3=derivative(Q2)
 minQ=min(v.lower() for v in bernstein(Q,a));maxQ=max(v.upper() for v in bernstein(Q,a))
 maxQ3=max(v.upper() for v in bernstein(Q3,a))
 check('actual_collar_elliptic',minQ>arb('.16') and maxQ<arb('.17'))
 check('actual_Q2_strictly_decreasing',maxQ3<arb('-3.8'))
 check('Q2_endpoint_signs',Q2[0]>0 and evaluate(Q2,a)<0)
 lo=arb('1e-12');hi=arb('2e-11')
 check('initial_root_bracket',evaluate(Q2,lo)>0 and evaluate(Q2,hi)<0)
 for iterations in range(180):
  mid=(lo+hi)/2;value=evaluate(Q2,mid)
  if value>0:lo=mid
  elif value<0:hi=mid
  else:break
 check('certified_root_bracket',evaluate(Q2,lo)>0 and evaluate(Q2,hi)<0)
 root=lo.union(hi)
 check('selected_root_interval',root>arb('1.0179e-11') and root<arb('1.0180e-11'))
 tinybern=bernstein(Q,hi);tinymin=min(v.lower() for v in tinybern);tinymax=max(v.upper() for v in tinybern)
 Y=arb(lo.lower()/tinymax.sqrt(),0).union(hi.upper()/tinymin.sqrt())
 check('proper_root_location',Y>arb('2.4933e-11') and Y<arb('2.4935e-11'))
 check('old_tiny_collar_not_extended',Y>arb('1e-15'))
 # Exact proper-coordinate chain rule for the actual auxiliary probe.
 z=S.symbols('z',positive=True);beta=S.symbols('beta',positive=True);q=S.Function('Q')(z)
 DY=lambda v:S.sqrt(q)*S.diff(v,z)
 vp=DY(S.sqrt(q))/S.sqrt(beta)
 check('exact_probe_identity',S.simplify(vp-S.diff(q,z)/(2*S.sqrt(beta)))==0)
 check('exact_probe_derivative',S.simplify(DY(vp)-S.diff(q,z,2)*S.sqrt(q)/(2*S.sqrt(beta)))==0)
 # Formal adjoint identity, derived by product rules on a physical core.
 A=S.Function('A')(z);B=S.Function('B')(z);V=S.Function('V')(z);f=S.Function('f')(z)
 adj=S.diff(A*f,z,2)-S.diff(B*f,z)+V*f
 target=A*S.diff(f,z,2)+(2*S.diff(A,z)-B)*S.diff(f,z)+(S.diff(A,z,2)-S.diff(B,z)+V)*f
 check('formal_adjoint_identity',S.expand(adj-target)==0)
 check('symmetric_divergence_identity',S.simplify((adj-(A*S.diff(f,z,2)+B*S.diff(f,z)+V*f)).subs(B,S.diff(A,z)).doit())==0)
 # Centered differences start with d^2 u'v'; their subtraction starts d^4.
 d=S.symbols('d');uj=S.symbols('u0:5');vj=S.symbols('v0:5')
 up=sum(uj[j]*(d/2)**j/S.factorial(j) for j in range(5));um=up.subs(d,-d)
 vv=sum(vj[j]*(d/2)**j/S.factorial(j) for j in range(5));vm=vv.subs(d,-d)
 rem=S.expand((up-um)*(vv-vm)-d*d*uj[1]*vj[1])
 check('symmetric_subtraction_fourth_order',all(rem.coeff(d,j)==0 for j in range(4)))
 check('symmetric_subtraction_coefficient',rem.coeff(d,4)==(uj[1]*vj[3]+uj[3]*vj[1])/24)
 # Exact strong row on z^s; differentiating the coefficients is mandatory.
 s=S.symbols('s');der=lambda expr:S.diff(expr,z)+s*expr/z
 local=der(A*s/z)+V/z
 row=S.expand(z*S.diff(z**4*local,z)+(s-4)*z**4*local)
 DA=z*S.diff(A,z);D2A=z*S.diff(DA,z);DV=z*S.diff(V,z)
 expected=z*z*(A*s*(s-1)*(s-2)+DA*(2*s*s-3*s)+D2A*s)+z**3*(V*(s-1)+DV)
 check('exact_weighted_divergence_row',S.expand(row-expected)==0)
 sig=S.Rational(801,200);t=S.symbols('t')
 polynomials=[s*(s-1)*(s-2),2*s*s-3*s,s]
 check('shifted_polynomials_positive',all(c>0 for poly in polynomials for c in S.Poly(poly.subs(s,sig+t),t).all_coeffs()))
 check('summed_contact_polynomial',S.expand(sum(polynomials)-s*s*(s-1))==0)
 check('summed_potential_polynomial',(s-1)+1==s)
 pp=S.symbols('p');hcut=[1,3,10,60]
 multiplier=sum(S.binomial(3,j)*sum(S.binomial(j,k)*hcut[k]*(pp+1)**(j-k) for k in range(j+1)) for j in range(4))
 M=lambda p:p**3+12*p*p+51*p+100
 check('coefficient_cutoff_product_constant',S.expand(multiplier-M(pp+1))==0)
 # For t=log(a/z)>=0, d[e^-pt(1+t)]/dt=e^-pt(1-p(1+t))<=0, p>=1.
 check('log_allowance_bound',S.simplify(S.diff(S.exp(-pp*t)*(1+t),t)-S.exp(-pp*t)*(1-pp-pp*t))==0)
 prior=json.loads((ROOT/'evidence/current-local-quantum-contact-graph.json').read_bytes())['runs'][-1]['result']
 caps=[old.decode(v) for v in prior['shape_multiplier_constants']];sigma=arb('4.005')
 TA=caps[3]+(3*sigma-1)*caps[2]+(3*sigma*sigma-2*sigma)*caps[1]+sigma*sigma*(sigma-1)*caps[0]
 TV=caps[1]+sigma*caps[0];fq=pars['C']/(arb.pi()*pars['B2'])
 fA=fq*a*a*388*TA;fV=fq*a**3*560*TV
 check('graph_sensitivity_guards',fA<arb('3.137593e-30') and fV<arb('3.746000e-33'))
 thresholdA=arb('.01')/fA;thresholdV=arb('.01')/fV
 check('contact_threshold_guards',thresholdA>arb('3.1871569908e27') and thresholdV>arb('2.6695143244e30'))
 selected=fA*arb('3.18e27')+fV*arb('2.66e30')
 check('conditional_contact_budget',selected<arb('.02'))
 check('conditional_combined_budget',selected+arb('.1599164')<arb('.1799164'))
 illustrative=fA*arb('1e27')+fV*arb('1e29')
 check('illustrative_budget',illustrative<arb('.003513'))
 return dict(Q_min=old.encoded(minQ),Q_max=old.encoded(maxQ),Q3_max=old.encoded(maxQ3),root_bracket=old.encoded(root),root_Y_enclosure=old.encoded(Y),root_iterations=iterations+1,
  root_low_Q2=old.encoded(evaluate(Q2,lo)),root_high_Q2=old.encoded(evaluate(Q2,hi)),two_C_root=old.encoded(2*pars['C']*root),
  actual_probe_division_uniformly_obstructed=True,old_tiny_collar_theorem_preserved=True,actual_operator_obstructed=False,
  symmetric_local_contact_identity_proved=True,full_symmetric_kernel_graph_bound_proved=False,
  fA=old.encoded(fA),fV=old.encoded(fV),TA=old.encoded(TA),TV=old.encoded(TV),physical_factor=old.encoded(fq),
  contact_threshold_A=old.encoded(thresholdA),contact_threshold_V=old.encoded(thresholdV),
  selected_MA='3.18e27',selected_MV='2.66e30',conditional_contact_upper=old.encoded(selected),illustrative_contact_upper=old.encoded(illustrative),
  proposed_nonlocal_target='.1599164',actual_contact_jets_certified=False,actual_nonlocal_remainder_certified=False,
  actual_data_target_certified=False,full_X_ball_admission_proved=False,full_signed_projection_executed=False,actual_controlled_solution=False,
  nonlinear_preserved_ball=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
  physical_source_or_norm_changed=False,normal_or_slow_root_changed=False,finite_physical_modes_selected=False)

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
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),summary={k:result[k] for k in ['root_bracket','root_Y_enclosure','fA','fV','conditional_contact_upper'] if k in result}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
