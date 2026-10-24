"""Actual-domain TT inverse-norm diagnostic; not the physical correction."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys, json, traceback
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,'C:/Users/dfred/AppData/Local/Temp/s169af-arb')
sys.path.insert(0,str(ROOT/'checks'))
import current_physical_york_dual_assembly as york
from flint import arb,ctx
STEM='current-weak-to-admission-obstruction'
OWN=['notes/'+STEM+'.md','checks/current_weak_to_admission_obstruction.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name, condition):
 ROWS.append({'name':name,'passed':bool(condition)})
 if not condition: raise AssertionError(name)
def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x): return x.str(35)
def fall(k,j): return factorial(k)//factorial(k-j) if k>=j else 0

class Majorant:
 """Absolute value, first and second derivative product majorants."""
 def __init__(self,v=0,d=0,dd=0): self.v,self.d,self.dd=map(arb,(v,d,dd))
 @staticmethod
 def lift(x):return x if isinstance(x,Majorant) else Majorant(x)
 def __add__(self,x):
  x=self.lift(x);return Majorant(self.v+x.v,self.d+x.d,self.dd+x.dd)
 __radd__=__add__
 def __mul__(self,x):
  x=self.lift(x);return Majorant(self.v*x.v,self.d*x.v+self.v*x.d,self.dd*x.v+2*self.d*x.d+self.v*x.dd)
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def __pow__(self,k):
  out=Majorant(1)
  for _ in range(k):out=out*self
  return out

def execute():
 p,h,zh,env=york.actual_parameters()
 qcoef={k:env[s] for k,s in zip(york.powers,york.qs)}
 qbound=[sum((abs(v).upper()*fall(k,j) for k,v in qcoef.items() if k>=j),arb(0)) for j in range(7)]
 z0=arb(9)/16;zc=arb(5)/8;length=arb(7)/16
 qcenter=[sum((v*fall(k,j)*zc**(k-j) for k,v in qcoef.items() if k>=j),arb(0)) for j in range(7)]
 check('support_before_actual_bolt',zh>arb(11)/16)
 check('actual_positive_center',qcenter[0]>0)
 n=arb('1e20');amp=n**(-arb(11)/2);jb=[amp*(n+832)**j for j in range(7)]
 # a=(Q'/z^4)' = Q''/z^4-4Q'/z^5, all its required jets.
 def powjet(power,j,z):
  value=arb(1)
  for k in range(j):value*=power-k
  return value*z**(power-j)
 ab=[];ac=[]
 for j in range(5):
  ab.append(sum((comb(j,k)*(qbound[k+2]*abs(powjet(-4,j-k,z0))+4*qbound[k+1]*abs(powjet(-5,j-k,z0))) for k in range(j+1)),arb(0)))
  ac.append(sum((comb(j,k)*(qcenter[k+2]*powjet(-4,j-k,zc)-4*qcenter[k+1]*powjet(-5,j-k,zc)) for k in range(j+1)),arb(0)))
 kb=[ab[0]*amp/8]
 for r in range(1,6):kb.append(sum((comb(r-1,j)*ab[j]*jb[r-1-j] for j in range(r)),arb(0)))
 nu=[];tt=[]
 for j in range(4):
  prod=sum((comb(j+1,k)*qbound[k]*jb[j+1-k] for k in range(j+2)),arb(0))
  tail=sum((comb(j,k)*fall(4,k)*kb[j-k] for k in range(min(j,4)+1)),arb(0))
  qjoverz=sum((comb(j,k)*abs(powjet(-1,k,z0))*sum((comb(j-k,l)*qbound[l]*jb[j-k-l] for l in range(j-k+1)),arb(0)) for k in range(j+1)),arb(0))
  nu.append(prod/2+tail/2+qjoverz);tt.append(3*prod/2+tail)
 C,B2,q=p['C'],p['B2'],arb(p['q']);la0=arb(2000)**-10
 Q,Q1,Q2=qbound[:3];iv=1/z0
 phi=Majorant(-z0.log(),nu[0]);phip=Majorant(iv,nu[1]);phipp=Majorant(iv**2,nu[2])
 LL=Majorant(Q2,2*(Q*tt[2]+Q1*tt[1]))
 MM=Q*phipp+Q1*phip;VV=Q*phip**2;DD=LL+6*MM+6*VV
 def expon(k,b):return Majorant(1,abs(k)*b,k*k*b*b)
 en2,en4,eu2,eum2=expon(2,nu[0]),expon(4,nu[0]),expon(2,tt[0]),expon(-2,tt[0])
 local={
  'Einstein':(en2*DD*iv**2/2+en2*eu2*iv**2/C+3*la0*B2/C*en4*eu2*iv**4).dd*length,
  'magnetic':(arb.pi()*q*q*arb(2000)**8/(C*B2)*eu2).dd*length,
  'R2':(8*arb.pi()*arb('1e22')/(C*B2)*(4*eu2+4*C*DD+C*C*eum2*DD**2)).dd*length,
  'Weyl':((C*C*eum2*(phi*LL**2/60+LL*MM/30+MM**2/10+LL*VV/30+7*MM*VV/90+7*VV**2/180)+C*phi*LL/18+(arb(1)/15+q*q/3)*phi*eu2+C*VV/18)/(2*arb.pi()*C*B2)).dd*length}
 U0=2*(604*q*q+(1+arb('6e-9'))*(q*q/6-arb(1)/30))/(arb.pi()*C*B2)
 delta=Q*tt[2]+Q1*tt[1]
 local['a0']=length*U0*tt[0]**2
 local['a2']=length*100*C/(480*arb.pi()*B2)*(4*Q2**2*tt[0]**2+16*Q2*tt[0]*delta+8*delta**2)
 local['LLL']=length*q/(6*arb.pi()*B2)*Q*tt[1]**2
 aux=arb(3)/2*Q1*tt[2]+Q2*tt[1]/2
 Z=length*Q*(arb(5)/8*Q**2*tt[3]**2+arb(3)/2*Q*tt[3]*aux+aux**2)
 s4=4*C*C/q
 local['complete_positive_geodesic']=2*s4*Z/(525*arb.pi()*B2)
 total=sum(local.values(),arb(0))/(1-arb('.0022'))
 weak=[]
 for N in [1,2]:
  W=N*arb.pi()*B2/p['gamma'];vnorm=(W*total).sqrt();dnorm=(arb('1.338')*W*total).sqrt()
  weak.append(dict(N=N,physical_action_abs=enc(W),V_upper=enc(vnorm),V_dual_residual_upper=enc(dnorm)))
  check('weak_diagnostic_under_1e_8_N'+str(N),dnorm<arb('1e-8'))
 # Differential operator coefficients, retaining all lower terms.
 qs=sp.symbols('q0:7',positive=True)
 coeff={1:sp.sqrt(qs[0])}
 for order in range(2,6):
  nxt={}
  for j,c in coeff.items():
   nxt[j]=nxt.get(j,0)+sp.sqrt(qs[0])*sum(sp.diff(c,qs[k])*qs[k+1] for k in range(6))
   nxt[j+1]=nxt.get(j+1,0)+sp.sqrt(qs[0])*c
  coeff={j:sp.expand(v) for j,v in nxt.items()}
 def evaluate(expr):
  if expr.is_Number:return arb(str(expr.p))/int(expr.q)
  if expr in qs:return qcenter[qs.index(expr)]
  if expr.is_Add:return sum((evaluate(x) for x in expr.args),arb(0))
  if expr.is_Mul:
   ans=arb(1)
   for x in expr.args:ans*=evaluate(x)
   return ans
  if expr.is_Pow:return evaluate(expr.base)**evaluate(expr.exp)
  raise ValueError(str(expr))
 jc=[amp*((-1)**(j//2))*n**j if j%2==0 else arb(0) for j in range(7)]
 kc=[arb(0,kb[0].upper())]
 for r in range(1,6):kc.append(sum((comb(r-1,k)*ac[k]*jc[r-1-k] for k in range(r)),arb(0)))
 tj=[]
 for j in range(6):
  # Remove the leading sixth J jet algebraically; independent interval
  # subtraction of two copies would lose the parameter correlation.
  prod=sum((comb(j+1,k)*qcenter[k]*jc[j+1-k] for k in range(1 if j==5 else 0,j+2)),arb(0))
  tail=sum((comb(j,k)*fall(4,k)*zc**(4-k)*kc[j-k] for k in range(min(j,4)+1)),arb(0))
  tj.append(-3*prod/2+tail)
 remainder=sum((evaluate(c)*tj[j] for j,c in coeff.items()),arb(0))
 leading=arb(3)/2*qcenter[0]**(arb(7)/2)*n.sqrt()
 fifth=leading+remainder
 check('fifth_derivative_positive',fifth>0)
 check('lower_terms_retained_and_small',abs(remainder)<leading*arb('1e-30'))
 diagnostic=fifth/(factorial(5)*100**5)
 check('proper_fifth_jet_diagnostic_above_1e_6',diagnostic>arb('1e-6'))
 bolt=(1+h)/(2*h)*zh**4*kb[0]
 check('bolt_response_below_physical_target',bolt<arb('1.667e-9'))
 # Independent exact finite algebra behind the cutoff and York constraint.
 t=sp.symbols('t');cut=sum(sp.binomial(13,k)*t**k*(1-t)**(13-k) for k in range(7))
 check('cutoff_beta_derivative',sp.expand(sp.diff(cut,t)+12012*t**6*(1-t)**6)==0)
 for j in range(1,7):check('cutoff_endpoint_jet_'+str(j),sp.diff(cut,t,j).subs(t,0)==0 and sp.diff(cut,t,j).subs(t,1)==0)
 z=sp.symbols('z',positive=True);Qf=sp.Function('Q')(z);J=sp.Function('J')(z);K=sp.Function('K')(z)
 xi=sp.diff(Qf*J,z)-z**4*K
 york_id=sp.diff(xi/z**4,z)-(Qf*sp.diff(J,z,2)+(2*sp.diff(Qf,z)-4*Qf/z)*sp.diff(J,z))/z**4
 check('global_actual_TT_identity',sp.simplify(york_id.subs(sp.diff(K,z),sp.diff(sp.diff(Qf,z)/z**4,z)*J))==0)
 return dict(n=enc(n),amplitude=enc(amp),C=enc(C),B2=enc(B2),h=enc(h),gamma=enc(p['gamma']),
  Q_center_jets=list(map(enc,qcenter)),Q_global_absolute_jet_bounds=list(map(enc,qbound)),
  a_global_jet_bounds=list(map(enc,ab)),K_global_upper=enc(kb[0]),
  canonical_nu_jet_bounds=list(map(enc,nu)),canonical_T_jet_bounds=list(map(enc,tt)),
  normalized_component_bounds={k:enc(v) for k,v in local.items()},weak=weak,
  complete_angular_fourth_inverse_mass_sum_upper=enc(s4),fifth_proper_derivative=enc(fifth),
  fifth_proper_leading_term=enc(leading),fifth_lower_jet_remainder=enc(remainder),unchanged_unit_fifth_jet_diagnostic=enc(diagnostic),
  fifth_operator_coefficients={str(k):str(v) for k,v in coeff.items()},physical_bolt_upper=enc(bolt),
  claimed_finite_V_dual_to_C5_constant=False,actual_correction_computed=False,
  nonlinear_admission_certified=False,physical_source_changed=False,all_sector_stability=False)

def main():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 paths=set(OWN+review['retained_inputs'])
 # Imported source modules are retained as inputs; no old checker is run.
 for module in list(sys.modules.values()):
  f=getattr(module,'__file__',None)
  if f:
   try:rel=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if rel.endswith('.py'):paths.add(rel)
 before={p:digest(p) for p in sorted(paths)}
 out=ROOT/'evidence'/f'{STEM}.json';archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',
  input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_text(json.dumps(archive,indent=2)+'\n',encoding='utf-8',newline='\n')
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
  run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if run['input_sha256_after']!=before:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 log=ROOT/'evidence'/f'{STEM}.log'
 with log.open('a',encoding='utf-8',newline='\n') as f:f.write(json.dumps(run)+'\n')
 print(json.dumps({k:v for k,v in run.items() if k not in ['input_sha256_before','input_sha256_after','frozen_sources_utf8']},indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
