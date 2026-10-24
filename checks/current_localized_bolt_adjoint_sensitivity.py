"""Full-domain scalar return and a cold-localized TT adjoint lower bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import factorial,comb
import json,sys,traceback
sys.path.insert(0,'C:/Users/dfred/AppData/Local/Temp/s169af-arb')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import current_weak_to_admission_obstruction as prev
from flint import arb,ctx
york=prev.york;Majorant=prev.Majorant
STEM='current-localized-bolt-adjoint-sensitivity'
OWN=['notes/'+STEM+'.md','checks/current_localized_bolt_adjoint_sensitivity.py','evidence/'+STEM+'-review.json']
ROWS=[];PARAM={};KP=None
def check(name,ok):
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def enc(x):return x.str(35)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def fall(k,j):return factorial(k)//factorial(k-j) if k>=j else 0

class R:
 """v^i z^j K^k with actual relation z=z_h-ha v and K'=KP."""
 def __init__(self,x=0):
  if isinstance(x,R):self.c=dict(x.c)
  elif isinstance(x,dict):self.c={k:arb(v) for k,v in x.items() if not arb(v).is_zero()}
  else:self.c={} if arb(x).is_zero() else {(0,0,0):arb(x)}
 def __add__(self,x):
  out=dict(self.c)
  for k,v in R(x).c.items():out[k]=out.get(k,arb(0))+v
  return R(out)
 __radd__=__add__
 def __neg__(self):return R({k:-v for k,v in self.c.items()})
 def __sub__(self,x):return self+-R(x)
 def __rsub__(self,x):return R(x)+-self
 def __mul__(self,x):
  out={}
  for k,c in self.c.items():
   for j,d in R(x).c.items():
    key=tuple(a+b for a,b in zip(k,j));out[key]=out.get(key,arb(0))+c*d
  return R(out)
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def __pow__(self,k):
  out=R(1)
  for _ in range(k):out=out*self
  return out
 def dv(self):
  out=R(0)
  for (i,j,k),c in self.c.items():
   if i:out+=R({(i-1,j,k):i*c})
   if j:out+=R({(i,j-1,k):-PARAM['ha']*j*c})
   if k:out+=R({(i,j,k-1):k*c})*KP
  return out
 def reduced(self):
  m=max(0,-min((k[1] for k in self.c),default=0));rows={}
  for (i,j,k),c in self.c.items():
   for p in range(j+m+1):
    key=(i+p,k)
    rows[key]=rows.get(key,arb(0))+c*comb(j+m,p)*PARAM['zh']**(j+m-p)*(-PARAM['ha'])**p
  return m,rows
 def bound(self,kmax,cells=8):
  m,rows=self.reduced();answer=arb(0)
  for power in sorted({k for i,k in rows}):
   coeff={i:c for (i,k),c in rows.items() if k==power};deg=max(coeff,default=0);best=arb(0)
   for cell in range(cells):
    lo=arb(cell)/cells;wid=arb(1)/cells
    transformed=[sum((c*comb(i,j)*lo**(i-j)*wid**j for i,c in coeff.items() if i>=j),arb(0)) for j in range(deg+1)]
    bern=[sum((transformed[j]*arb(comb(k,j))/comb(deg,j) for j in range(k+1)),arb(0)) for k in range(deg+1)]
    for b in bern:
     if b.abs_upper()>best:best=b.abs_upper()
   answer+=best*kmax**power
  return answer/PARAM['za']**m
 def serial(self):return [{'v':i,'z':j,'K':k,'coefficient':enc(c)} for (i,j,k),c in sorted(self.c.items())]

def mon(i=0,j=0,k=0):return R({(i,j,k):1})
class Dual:
 def __init__(self,v=0,d=0):self.v,self.d=R(v),R(d)
 @staticmethod
 def lift(x):return x if isinstance(x,Dual) else Dual(x)
 def __add__(self,x):x=self.lift(x);return Dual(self.v+x.v,self.d+x.d)
 __radd__=__add__
 def __neg__(self):return Dual(-self.v,-self.d)
 def __sub__(self,x):return self+-self.lift(x)
 def __rsub__(self,x):return self.lift(x)+-self
 def __mul__(self,x):x=self.lift(x);return Dual(self.v*x.v,self.d*x.v+self.v*x.d)
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def __pow__(self,k):
  out=Dual(1)
  for _ in range(k):out=out*self
  return out
 def dv(self):return Dual(self.v.dv(),self.d.dv())

def execute():
 global KP
 p,h,zh,env=york.actual_parameters();a=arb('1e-4');ha=h*a;za=zh-ha
 PARAM.update(h=h,a=a,ha=ha,zh=zh,za=za)
 C,B2,lam,qmag=p['C'],p['B2'],p['g'],arb(p['q']);beta=lam*B2/C
 check('actual_negative_d',p['eps']<0)
 check('global_P_lower_and_monotonicity',p['beta']>arb('.16') and 3*p['beta']+p['u']+p['eps']>0)
 coeff={k:env[s] for k,s in zip(york.powers,york.qs)}
 qc=[arb(0),2*p['gamma']/h]
 for j in range(2,7):qc.append((-1)**j*h**(j-2)*sum((c*fall(k,j)*zh**(k-j) for k,c in coeff.items() if k>=j),arb(0))/factorial(j))
 q=sum((mon(i=j)*c*a**j for j,c in enumerate(qc)),R(0));psi=(1-mon(i=1))**8
 Apr=(q.dv().dv()/a**2)*mon(j=-4)+4*h*(q.dv()/a)*mon(j=-5)
 KP=-a*Apr*psi
 Amax=Apr.bound(arb(0));Kmax=a*Amax/9
 check('A_positive_on_support',qc[2]>arb('.9') and q.bound(arb(0))<arb('.00021'))
 # q'' positive and q' positive follow from positive low coefficients
 # and the absolute high-degree coefficient budget.
 high2=sum((abs(qc[j])*j*(j-1)*a**(j-2) for j in range(3,7)),arb(0))
 high1=sum((abs(qc[j])*j*a**(j-1) for j in range(2,7)),arb(0))
 check('q_second_positive',2*qc[2]-high2>arb('1.9'))
 check('q_first_positive',qc[1]-high1>arb('1.9'))
 xi=-(q*psi).dv()/a-mon(j=4,k=1)
 nu=xi/2+h*q*psi*mon(j=-1);T=-xi+(q*psi).dv()/(2*a)
 # Exact full local scalar/mixed density, in the regular scaled variable.
 phiv=Dual(ha*mon(j=-1),nu.dv());ushape=Dual(0,T)
 lap=lambda f:(f.dv()*q).dv()/a**2
 LL=Dual(-q.dv().dv()/a**2)-2*lap(ushape)
 MM=(phiv*q).dv()/a**2;VV=phiv**2*q/a**2;DD=LL-6*MM-6*VV
 en2,en4,eu2,em2=Dual(1,2*nu),Dual(1,4*nu),Dual(1,2*T),Dual(1,-2*T)
 UU=2+C*em2*DD;la0=arb(2000)**-10
 FE=en2*DD*mon(j=-2)+2*en2*eu2*mon(j=-2)/C+12*la0*B2/C*en4*eu2*mon(j=-4)
 FR=-96*arb.pi()*arb('1e22')/B2*(lap(UU)-2*(UU*phiv*q).dv()/a**2)
 FW=(C*C*(em2*LL**2/60+lap(em2*LL)/30+lap(em2*(-MM/5-7*VV/90))-2*(em2*(LL/30-7*MM/90-7*VV/90)*phiv*q).dv()/a**2)-C*LL/18+(arb(1)/15-qmag*qmag/3)*eu2+C*MM/9)/(2*arb.pi()*C*B2)
 g=-(FE.d+FR.d+FW.d)*mon(j=4)/beta
 check('mixed_source_linear_in_global_K',all(k<=1 for i,j,k in g.c))
 gmax=g.bound(Kmax,cells=16)
 kappa=2/(1+h);I=(1+a/kappa).log()+(a/kappa)*(1+kappa/a).log()
 VP=zh**2/(p['beta']*za**4)*I;m=1-arb('1.1e-8');eps=arb('1e-26')
 us=gmax*(beta*VP/(6*m)+eps/(6*arb('.63')*m*(1-eps)))
 blow=p['gamma']/h-us
 check('complete_mixed_source_bound',gmax<100)
 check('actual_scalar_bolt_return_small',us<arb('.02'))
 check('effective_bolt_functional_nonzero',blow>arb('.98'))
 # All positive local Hessian majorants use t derivatives, NOT z jets.
 nj=[];tj=[];qj=[];nv,ttv,qv=nu,T,q
 for j in range(4):
  nj.append(nv.bound(Kmax)/a**j);tj.append(ttv.bound(Kmax)/a**j);qj.append(qv.bound(arb(0))/a**j)
  nv,ttv,qv=nv.dv(),ttv.dv(),qv.dv()
 Q,Q1,Q2=qj[:3];iv=1/za;length=ha
 phi=Majorant(-za.log(),nj[0]);phit=Majorant(h*iv,nj[1]);phitt=Majorant(h*h*iv**2,nj[2])
 L=Majorant(Q2,2*(Q*tj[2]+Q1*tj[1]));M=Q*phitt+Q1*phit;V=Q*phit**2;D=L+6*M+6*V
 ex=lambda k,b:Majorant(1,abs(k)*b,k*k*b*b)
 e2,e4,b2,bm2=ex(2,nj[0]),ex(4,nj[0]),ex(2,tj[0]),ex(-2,tj[0])
 local={
  'Einstein':length*(e2*D*iv**2/2+e2*b2*iv**2/C+3*la0*B2/C*e4*b2*iv**4).dd,
  'magnetic':length*(arb.pi()*qmag*qmag*arb(2000)**8/(C*B2)*b2).dd,
  'R2':length*(8*arb.pi()*arb('1e22')/(C*B2)*(4*b2+4*C*D+C*C*bm2*D**2)).dd,
  'Weyl':length*((C*C*bm2*(phi*L**2/60+L*M/30+M**2/10+L*V/30+7*M*V/90+7*V**2/180)+C*phi*L/18+(arb(1)/15+qmag*qmag/3)*phi*b2+C*V/18)/(2*arb.pi()*C*B2)).dd}
 U0=2*(604*qmag*qmag+(1+arb('6e-9'))*(qmag*qmag/6-arb(1)/30))/(arb.pi()*C*B2)
 delta=Q*tj[2]+Q1*tj[1]
 local['a0']=length*U0*tj[0]**2
 local['a2']=length*100*C/(480*arb.pi()*B2)*(4*Q2**2*tj[0]**2+16*Q2*tj[0]*delta+8*delta**2)
 local['LLL']=length*qmag/(6*arb.pi()*B2)*Q*tj[1]**2
 aux=arb(3)/2*Q1*tj[2]+Q2*tj[1]/2
 Z=length*Q*(arb(5)/8*Q**2*tj[3]**2+arb(3)/2*Q*tj[3]*aux+aux**2)
 s4=4*C*C/qmag;local['full_positive_geodesic']=2*s4*Z/(525*arb.pi()*B2)
 eta=arb('.0022');MS=arb('1.341');ell=(1+h)/h;Lc=sum(local.values(),arb(0))/(1-eta)
 norms=[]
 for N in [1,2]:
  W=N*arb.pi()*B2/p['gamma'];dual=ell*blow/(W*Lc).sqrt();adjoint=dual/MS
  norms.append(dict(N=N,W_abs=enc(W),effective_functional_dual_lower=enc(dual),comparison_adjoint_norm_lower=enc(adjoint)))
  check('actual_adjoint_norm_large_N'+str(N),adjoint>arb('1e8'))
 old=json.loads((ROOT/'evidence/current-schur-bolt-functional.json').read_bytes())['runs'][-1]
 glob=next(x for x in old['numerics'] if x['kind']=='actual_Schur_trial_diagonal_and_source_threshold' and x['name']=='TT_zero_J')
 Lg=york.endpoint.raw.decode(glob['actual_l_diagonal_upper_over_action'])
 threshold=arb('1.667e-9')*MS**2*(Lg*Lc).sqrt()/(eta*(1+eta)*ell*blow)
 check('new_necessary_screen_below_1e_52',threshold<arb('1e-52') and threshold>0)
 allowance=arb('1.3079691e-40');rest=arb('2.123e-42')
 rhs=arb('1.667e-9')*MS**2*(Lg*Lc).sqrt()/(ell*blow*allowance)
 hypothetical_eta=2*rhs/(1+(1+4*rhs).sqrt())
 gm,grows=g.reduced()
 return dict(parameters={k:enc(v) for k,v in dict(h=h,zh=zh,gamma=p['gamma'],a=a,C=C,B2=B2,beta=beta).items()},
  exact_scaled_q_t_coefficients=list(map(enc,qc)),A_upper=enc(Amax),K_upper=enc(Kmax),
  complete_mixed_source_representation=g.serial(),source_denominator_power=gm,
  reduced_source_coefficients=[dict(v=i,K=k,coefficient=enc(c)) for (i,k),c in sorted(grows.items())],
  g_max=enc(gmax),poisson_factor=enc(VP),scalar_return_upper=enc(us),effective_point_lower=enc(blow),
  nu_t_jet_bounds=list(map(enc,nj)),T_t_jet_bounds=list(map(enc,tj)),q_t_jet_bounds=list(map(enc,qj)),
  normalized_energy_component_bounds={k:enc(v) for k,v in local.items()},cold_energy_upper_over_action=enc(Lc),
  broad_trial_energy_upper_over_action=enc(Lg),norms=norms,necessary_projection_ceiling=enc(threshold),
  existing_full_source_allowance_to_ceiling=enc(allowance/threshold),existing_rest_allowance_to_ceiling=enc(rest/threshold),
  hypothetical_eta_if_full_allowance_were_actual_size=enc(hypothetical_eta),
  actual_source_sign_measured=False,actual_source_lower_bound_proved=False,uniform_certificate_ruled_out_for_actual_source=False,
  actual_correction_computed=False,nonlinear_admission=False,nonexistence=False,all_sector_stability=False)

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
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
  run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{});print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),
  summary={k:result.get(k) for k in ['g_max','scalar_return_upper','effective_point_lower','cold_energy_upper_over_action','norms','necessary_projection_ceiling','existing_full_source_allowance_to_ceiling','hypothetical_eta_if_full_allowance_were_actual_size']}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
