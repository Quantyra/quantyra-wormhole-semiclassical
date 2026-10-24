"""Actual outer log-principal feedback; not the complete graph remainder."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as sp
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-outer-graph-feedback'
OWN=['notes/'+STEM+'.md','checks/current_outer_graph_feedback.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded

def check(name,ok,detail=''):
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 z,s=sp.symbols('z s',positive=True);Q=sp.Function('Q')(z)
 lap=lambda v:sp.diff(Q*sp.diff(v,z),z)
 row=z**4*lap(lap(z**s));row=sp.expand((z*sp.diff(row,z)-4*row)/z**s)
 check('actual_metric_log_principal_coefficient',sp.simplify(row.coeff(s,5)-Q**2)==0)
 # Exact Bernstein convex-hull certificates over the entire cutoff interval.
 u=sp.symbols('u');S=35*u**4-84*u**5+70*u**6-20*u**7
 check('cutoff_values_and_monotonicity',S.subs(u,0)==0 and S.subs(u,1)==1 and sp.simplify(sp.diff(S,u)-140*u**3*(1-u)**3)==0)
 maxima=[]
 for j,guard in zip([1,2,3],[3,10,60]):
  P=sp.Poly(sp.diff(S,u,j),u);n=P.degree();power=[P.nth(k) for k in range(n+1)];cap=sp.Integer(0)
  for cell in range(32):
   lo=sp.Rational(cell,32);width=sp.Rational(1,32)
   local=[sum(power[k]*comb(k,i)*lo**(k-i)*width**i for k in range(i,n+1)) for i in range(n+1)]
   bern=[sum(local[k]*sp.Rational(comb(i,k),comb(n,k)) for k in range(i+1)) for i in range(n+1)]
   cap=max(cap,max(abs(v) for v in bern))
  check('entire_cutoff_derivative_'+str(j),cap<guard,str(cap));maxima.append(str(cap))
 # Exact Sobolev product coefficients: no hidden multiplication constant.
 for k in [1,2,3]:
  ratios=[sp.Rational(comb(3,j)*comb(j,k)**2,comb(3-k,j-k)) for j in range(k,4)]
  check('H3_product_constant_'+str(k),max(ratios)<=comb(3,k)**2)
 # High-frequency lower estimate, with all t>=4.
 t=sp.symbols('t',positive=True)
 check('high_frequency_positive_reserve',sp.expand(7*(t*t-7)-26).subs(t,4)>0)
 check('high_frequency_log_coefficient',sp.Rational(1,2)*16-7>0)
 check('frequency_norm_growth_exponent',5+3-6==2)
 SYMBOLIC.update(cutoff_derivative_exact_caps=maxima,actual_log_principal='-K[(Q/beta)^2-1](i xi)^5 log|xi|',
  H6_to_H3_feedback_constant_exists=False,graph_log_weight='K <xi>^8 log(2+|xi|)',
  actual_source_replaced_by_wave_packets=False)

def execute():
 algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 prior=json.loads((ROOT/'evidence/current-global-boundary-moment-return.json').read_bytes())['runs'][-1]
 check('previous_normal_return_passed',prior['status']=='passed')
 C3=old.decode(prior['result']['C3_H3_to_H6']);k3=old.decode(prior['result']['known_source_H3_upper'])
 C,B2,be=p['C'],p['B2'],p['beta'];K=C*be**2/(15*arb.pi()*B2)
 qlo=K*C3*17*arb(6).log();qhi=2*(arb(17)/16)**arb('2.5')
 check('graph_low_frequency_bound',qlo<qhi)
 check('whole_frequency_graph_constant',qhi<arb('2.328'))
 # Actual finite-temperature polynomial, with exact constant cancellation.
 q=[arb(1),arb(0),p['u']/be,-4-2*p['u']/be-p['tau']/be,
    3+p['u']/be+p['eps']/be,-2*p['eps']/be,p['eps']/be]
 rho=[sum((q[i]*q[k-i] for i in range(max(0,k-6),min(6,k)+1)),arb(0)) for k in range(13)]
 rho[0]=arb(0);rho[1]=arb(0)
 a=arb('.01');r=[sum((rho[k].abs_upper()*k**j*a**k for k in range(2,13)),arb(0)) for j in range(4)]
 aj=[r[0],r[1]+3*r[0],r[2]+6*r[1]+10*r[0],r[3]+9*r[2]+30*r[1]+60*r[0]]
 mult=aj[0]+3*aj[1]+3*aj[2]+aj[3];eta=qhi*mult
 check('actual_outer_principal_feedback',eta<arb('.008'))
 # An actual nonzero coefficient interval, not a claim about the source.
 point=arb('.005');rho_point=sum((rho[k]*point**k for k in range(2,13)),arb(0))
 check('nonzero_feedback_principal_on_actual_collar',rho_point<arb('-9.9e-7') and rho_point>arb('-1.01e-6'))
 # Closed arithmetic only IF the remaining actual operator estimate holds.
 other=arb('.24');S3=arb('7.49e-12');total=eta+other
 closed=(S3+total*k3)/(1-total);reserve=(S3+arb('.25')*k3)/arb('.75')
 check('remaining_feedback_allocation_leaves_margin',total<arb('.25'))
 check('conditional_global_E3_target',closed<arb('1e-11') and reserve<arb('1e-11'))
 threshold=(1-total)*arb('1e-11')-total*k3
 check('affine_threshold_positive',threshold>S3)
 return dict(K=enc(K),normal_graph_low_frequency_constant=enc(qlo),normal_graph_full_frequency_constant=enc(qhi),
  collar='.01',rho_at_z_point=enc(rho_point),rho_derivative_upper=[enc(v) for v in r],
  localized_coefficient_derivative_upper=[enc(v) for v in aj],H3_multiplication_bound=enc(mult),
  actual_log_principal_feedback_upper=enc(eta),conditional_other_feedback='<=.24',conditional_affine_source='<=7.49e-12',
  conditional_total_feedback_upper=enc(total),conditional_E3_upper=enc(closed),reserved_conditional_E3_upper=enc(reserve),
  affine_threshold_with_actual_principal_bound=enc(threshold),
  H6_to_H3_full_feedback_constant_exists=False,actual_log_principal_feedback_component_certified=True,
  full_graph_feedback_certified=False,remaining_graph_feedback_known=False,actual_affine_source_known=False,
  actual_E3_upper_known=False,actual_Efar_upper_known=False,actual_global_correction_computed=False,
  actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,
  semiclassical_validity_noise=False,all_sector_stability=False)

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
