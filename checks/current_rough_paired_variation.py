"""Exact full-frequency tangent test; NOT a curved-family heat theorem.

Run after reading notes/current-rough-paired-variation.md and its review.
All executions, source snapshots, input hashes and failures are retained.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json
import platform
import sys
import traceback

import sympy as sp
sys.path.insert(0, str(Path.home()/'AppData/Local/Temp/s169af-arb'))
import flint
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-rough-paired-variation'
INPUTS=[
 'checks/current_rough_paired_variation.py',
 'notes/current-rough-paired-variation.md',
 'evidence/current-rough-paired-variation-review.json',
 'notes/current-signed-resolvent-projection.md',
 'evidence/current-signed-resolvent-projection.json',
 'evidence/current-signed-resolvent-projection-closeout-verification.json',
 'notes/current-paired-regulator-domain-audit.md',
 'notes/current-smeared-heat-projection-domain.md',
 'notes/stationary-paired-dirac-phase.md',
 'notes/stationary-massive-boundary-layer.md',
 'notes/finite-outer-admission-and-response.md',
 'notes/finite-cold-admission-and-response.md',
 'notes/finite-comparison-conormal-regularity.md',
 'notes/current-controlled-solution-audit-closeout.md',
]

def hashes():
 return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}

def execute(check, emit):
 ctx.prec=512
 a,x,z,y,s=sp.symbols('a x z y s', positive=True)
 t=sp.symbols('t', real=True)
 aa=a*(1-a)
 icoeff=[]
 for n in range(9):
  pol=(-aa)**n/sp.factorial(n)
  if n: pol+=(1-2*a)**2/4*(-aa)**(n-1)/sp.factorial(n-1)
  value=sp.integrate(pol,(a,0,1))
  expected=1 if n==0 else (-1)**n*sp.factorial(n)**2/(2*sp.factorial(n)*sp.factorial(2*n+1))
  check('complete_alpha_coefficient_'+str(n),sp.simplify(value-expected)==0)
  icoeff.append(value)
 fc=[(1-n)*v for n,v in enumerate(icoeff)]
 check('F_local_coefficients',fc[:5]==[1,0,-sp.Rational(1,120),sp.Rational(1,840),-sp.Rational(1,10080)])
 # Elementary identity supplies a pointwise bound on ALL nonnegative x.
 phi=(1+z)*sp.exp(-z)
 rem=phi-1+z*z/2
 check('positive_remainder_derivative',sp.simplify(sp.diff(rem,z)-z*(1-sp.exp(-z)))==0)
 check('cubic_upper_moment',sp.integrate(aa**3/6,(a,0,1))==sp.Rational(1,840))
 check('quadratic_upper_moment',sp.integrate(aa**2/4,(a,0,1))==sp.Rational(1,120))
 # Four ordinary integrations by parts; h=O(z^4) at zero and O(z^3)
 # at infinity. Each boundary term vanishes. No finite part is used.
 h=phi-1+z*z/2-z**3/3
 check('IBP_fourth_derivative',sp.simplify(sp.diff(h,z,4)-(z-3)*sp.exp(-z))==0)
 check('IBP_origin_fourth_order',all(sp.limit(sp.diff(h,z,j),z,0)==0 for j in range(4)))
 weight=sp.prod(sp.Rational(2,2*j+1) for j in range(4))
 K=sp.simplify(weight*(sp.gamma(sp.Rational(3,2))-3*sp.gamma(sp.Rational(1,2))))
 beta=sp.simplify(sp.gamma(sp.Rational(9,2))**2/sp.gamma(9))
 J=sp.simplify(K*beta/4)
 check('J_exact',sp.simplify(J+5*sp.pi**sp.Rational(3,2)/49152)==0)
 cubic=sp.simplify(144*J/sp.pi**2)
 mixed=sp.simplify((-288*sp.pi)*fc[4]*(-sp.Rational(1,2))/sp.pi**2)
 check('cubic_coefficient',sp.simplify(cubic+15/(1024*sp.sqrt(sp.pi)))==0)
 check('mixed_log_coefficient',mixed==-1/(70*sp.pi))
 fourier=-2*sp.gamma(s+1)*sp.sin(sp.pi*s/2)
 check('mass_fourier_tail',fourier.subs(s,3)==12)
 check('log_fourier_tail',sp.simplify(sp.diff(fourier,s).subs(s,4)+24*sp.pi)==0)
 # Half-line Parseval: half the full-line trace, two Fourier signs,
 # and full-line H3 norm = sqrt(2) times the physical-half norm.
 check('physical_half_H3_factor',sp.Rational(1,2)/sp.pi**2*sp.Rational(1,840)*(2*sp.pi)*2==1/(420*sp.pi))
 # Exact collar Gram matrix includes the mass/log cross term.
 c,d,delta=sp.symbols('c d delta', real=True, positive=True)
 jet=sp.diff(c*y**3+d*y**4*sp.log(y),y,3)
 check('third_jet',sp.expand(jet-6*c-d*y*(24*sp.log(y)+26))==0)
 gram=36*delta*c*c+delta**2*(144*sp.log(delta)+84)*c*d+delta**3*(192*sp.log(delta)**2+288*sp.log(delta)+sp.Rational(388,3))*d*d
 check('collar_Gram_antiderivative',sp.simplify(sp.diff(gram,delta)-jet.subs(y,delta)**2)==0)
 check('collar_Gram_lower_endpoint',sp.limit(gram,delta,0)==0)
 # Uniform log regularization polynomial bounds over t in [0,1].
 E=[0,t]
 for j in range(1,5): E.append(sp.expand(2*t*(1-t)*sp.diff(E[-1],t)-j*E[-1]))
 fall=lambda n,k:factorial(n)//factorial(n-k) if k<=n else 0
 guards=[]
 for j in range(1,6):
  P=sp.expand(sum(comb(j,k)*fall(4,j-k)*E[k] for k in range(1,j+1)))
  guard=sum(abs(q) for q in sp.Poly(P,t).all_coeffs())+fall(4,j)
  check('weighted_log_guard_'+str(j),guard<=factorial(j)*100**j)
  guards.append({'j':j,'E':str(E[j]),'P':str(P),'coefficient_l1_guard':str(guard),'allowed':factorial(j)*100**j})
 # Elementary exact moments: diagnostics, not a proof of the curved bound.
 check('bridge_beta',sp.integrate(a*a*(1-a)**2,(a,0,1))==sp.Rational(1,30))
 spatial=delta**3*((1-sp.log(delta))/3+sp.Rational(1,9))
 check('dyadic_weight_moment',sp.simplify(sp.diff(spatial,delta)-delta**2*(1-sp.log(delta)))==0)
 pi=arb.pi(); B=arb('.001'); r=arb('.125'); loss=(2*B).exp()*(r.exp()-1)
 resolvent=(2*B+r).exp()/(1-loss)
 check('complex_neumann_margin',loss<arb('.134'))
 check('complex_resolvent_bound',resolvent<arb('1.32'))
 check('first_derivative_bound',resolvent/r<arb('10.56'))
 check('second_derivative_bound',4*resolvent/r**2<arb('337.92'))
 deltaa=arb('.1'); logd=deltaa.log(); t0=arb('.001')
 G00=36*deltaa;G01=deltaa**2*(72*logd+42)
 G11=deltaa**3*(192*logd**2+288*logd+arb(388)/3)
 det=G00*G11-G01**2
 check('collar_Gram_positive',G00>0 and det>0)
 check('curved_Dini_factor',2*t0.sqrt()<arb('.063246'))
 check('retained_sufficient_screen',arb('2.123e-42')<arb('1.287e-40'))
 check('retained_source_sum',sum(arb(v) for v in ['1.284506e-40','2e-43','2.123e-42','2.331e-44'])<arb('2.60e-40'))
 emit('exact',{'F_coefficients':[str(v) for v in fc], 'J':str(J),'cubic':str(cubic),'mixed_log':str(mixed),'H3_factor':'1/(420*pi)','regularization_guards':guards})
 emit('intervals', {k:str(v) for k,v in {
  'cubic_coefficient':-15/(1024*pi.sqrt()),'mixed_log_coefficient':-1/(70*pi),
  'tangent_integrated_H3_factor':t0/(420*pi),'curved_Dini_factor':2*t0.sqrt(),
  'collar_Gram_00':G00,'collar_Gram_01':G01,'collar_Gram_11':G11,'collar_Gram_determinant':det,
  'dyadic_weight_moment':deltaa**3*((1-logd)/3+arb(1)/9),
  'complex_neumann_loss':loss,'complex_resolvent_prefactor':resolvent,
  'first_metric_variation_prefactor':resolvent/r,'second_metric_variation_prefactor':4*resolvent/r**2,
  'Km_threshold':arb('1e-9')/arb('1.632e-11'),'Kc_threshold':arb('2e-9')/arb('1.36e-9'),
  'G_threshold_if_ED_point025':arb(800)/201,
 }.items()})

def main():
 output=ROOT/'evidence'/f'{STEM}.json'; log=ROOT/'evidence'/f'{STEM}.log'
 archive=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
 run={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running',
  'dependencies':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':512},
  'input_sha256_before':hashes(),'checks':[],'numerics':[],
  'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in INPUTS[:3]},
  'rough_curved_C2_action_proved':False,'actual_controlled_solution':False,
  'full_goal':'stable semiclassical wormhole; active and incomplete'}
 archive['runs'].append(run)
 def save():output.write_text(json.dumps(archive,indent=2)+'\n',encoding='utf-8',newline='\n')
 def row(obj):
  with log.open('a',encoding='utf-8',newline='\n') as f:f.write(json.dumps(obj)+'\n')
 def check(label,condition):
  item={'label':label,'passed':bool(condition)};run['checks'].append(item);row(item)
  if not condition:raise AssertionError(label)
 def emit(kind,value):run['numerics'].append({'kind':kind,**value});row({'kind':kind,**value})
 save();row({'run':len(archive['runs']),'started_utc':run['started_utc']})
 try:
  execute(check,emit)
  check('input_identity_preserved',hashes()==run['input_sha256_before'])
  run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc();row({'exception':run['exception']})
 finally:
  run['input_sha256_after']=hashes();run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 print(json.dumps({'status':run['status'],'checks':len(run['checks']),'inputs':len(INPUTS),'runs_preserved':len(archive['runs']),'rough_curved_C2_action_proved':False,'actual_controlled_solution':False}))
 if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__': main()
