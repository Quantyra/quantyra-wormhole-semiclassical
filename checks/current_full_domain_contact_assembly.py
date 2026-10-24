"""Whole-domain local-contact assembly; complete low-angular trace stays open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import json,sys,traceback,heapq
sys.path.insert(0,str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import sympy as sp
import current_positive_time_resolvent_reduction as prior
root=Path(__file__).resolve().parents[1];gy=prior.gy;bh=prior.bh
encoded=prior.encoded;decode=prior.previous.decode
STEM='current-full-domain-contact-assembly';PREVIOUS=prior.STEM
OWN=['notes/'+STEM+'.md','checks/current_full_domain_contact_assembly.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];POWER=arb(0)
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS},name
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name+': '+str(detail))
gy.check=check;bh.prior.scalar.prev.prior.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def hull(lo,hi):
 lo=arb(lo.lower());hi=arb(hi.upper());return (lo+hi)/2+arb(0,(hi-lo)/2)

class FP:
 """Finite z^k y^(j+t POWER) log(z)^l sum, y=1-z."""
 def __init__(self,x=0):
  if isinstance(x,FP):self.c=dict(x.c)
  elif isinstance(x,dict):self.c={k:arb(v) for k,v in x.items() if not arb(v).is_zero()}
  else:self.c={} if arb(x).is_zero() else {(0,0,0,0):arb(x)}
 def __add__(self,x):
  if isinstance(x,(LF,Jet)):return NotImplemented
  x=FP(x);d=dict(self.c)
  for k,v in x.c.items():d[k]=d.get(k,arb(0))+v
  return FP(d)
 __radd__=__add__
 def __neg__(self):return FP({k:-v for k,v in self.c.items()})
 def __sub__(self,x):
  if isinstance(x,(LF,Jet)):return NotImplemented
  return self+-FP(x)
 def __rsub__(self,x):return FP(x)+-self
 def __mul__(self,x):
  if isinstance(x,(LF,Jet)):return NotImplemented
  x=FP(x);d={}
  for k,v in self.c.items():
   for j,w in x.c.items():
    key=tuple(a+b for a,b in zip(k,j));d[key]=d.get(key,arb(0))+v*w
  return FP(d)
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def __pow__(self,n):
  assert isinstance(n,int) and n>=0
  out=FP(1)
  for _ in range(n):out=out*self
  return out
 def diff(self):
  out={}
  def add(k,c):out[k]=out.get(k,arb(0))+c
  for (k,j,t,l),c in self.c.items():
   if k:add((k-1,j,t,l),k*c)
   if j or t:add((k,j-1,t,l),-(j+t*POWER)*c)
   if l:add((k-1,j,t,l-1),l*c)
  return FP(out)
 def value(self,z):
  y=1-z;lg=z.log();yp=y**POWER;zp={};yq={};out=arb(0)
  for (k,j,t,l),c in self.c.items():
   if k not in zp:zp[k]=z**k
   if (j,t) not in yq:yq[j,t]=y**j*yp**t
   out+=c*zp[k]*yq[j,t]*lg**l
  return out
 def absrange(self,lo,hi):
  out=arb(0);zp={};yq={};lg=abs(lo.log())
  for (k,j,t,l),c in self.c.items():
   if k not in zp:zp[k]=arb(max((lo**k).upper(),(hi**k).upper()))
   if (j,t) not in yq:
    e=j+t*POWER;yq[j,t]=arb(max(((1-hi)**e).upper(),((1-lo)**e).upper()))
   out+=arb(abs(c).upper())*zp[k]*yq[j,t]*lg**l
  return out
 def integral(self,y0):
  out=arb(0)
  for (k,j,t,l),c in self.c.items():
   assert k>=0 and l==0,(k,j,t,l)
   m=arb(0)
   for i in range(k+1):
    a=j+t*POWER+i+1
    v=-y0.log() if a.is_zero() else -(a*y0.log()).expm1()/a
    m+=(-1)**i*comb(k,i)*v
   out+=c*m
  return out
 def serial(self):return [dict(powers=list(k),coefficient=encoded(v)) for k,v in sorted(self.c.items())]

def mon(k=0,j=0,t=0,l=0):return FP({(k,j,t,l):1})

class LF:
 """A linear differential expression in n,u, plus a known term."""
 def __init__(self,x=0):
  if isinstance(x,LF):self.c={k:FP(v) for k,v in x.c.items()}
  elif isinstance(x,dict):self.c={k:FP(v) for k,v in x.items() if FP(v).c}
  else:self.c={} if not FP(x).c else {(None,0):FP(x)}
 def __add__(self,x):
  if isinstance(x,Jet):return NotImplemented
  x=LF(x);d={k:FP(v) for k,v in self.c.items()}
  for k,v in x.c.items():d[k]=d.get(k,FP())+v
  return LF(d)
 __radd__=__add__
 def __neg__(self):return LF({k:-v for k,v in self.c.items()})
 def __sub__(self,x):return self+-LF(x)
 def __rsub__(self,x):return LF(x)+-self
 def __mul__(self,x):
  if isinstance(x,Jet):return NotImplemented
  assert not isinstance(x,LF),'A nonlinear field expression is not permitted'
  return LF({k:v*FP(x) for k,v in self.c.items()})
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def diff(self):
  out=LF()
  for (field,order),c in self.c.items():
   out+=LF({(field,order):c.diff()})
   if field is not None:out+=LF({(field,order+1):c})
  return out
 def strong(self):
  out={'n':FP(),'u':FP(),'known':FP()}
  for (field,order),c in self.c.items():
   if field is None:out['known']+=c;continue
   for _ in range(order):c=-c.diff()
   out[field]+=c
  return out

def atom(name):return LF({(name,0):FP(1)})

class Jet:
 def __init__(self,v=0,i=0,j=0,ij=0):self.v,self.i,self.j,self.ij=FP(v),LF(i),FP(j),LF(ij)
 @staticmethod
 def lift(x):return x if isinstance(x,Jet) else Jet(x)
 def __add__(self,x):
  x=Jet.lift(x);return Jet(self.v+x.v,self.i+x.i,self.j+x.j,self.ij+x.ij)
 __radd__=__add__
 def __neg__(self):return Jet(-self.v,-self.i,-self.j,-self.ij)
 def __sub__(self,x):return self+-Jet.lift(x)
 def __rsub__(self,x):return Jet.lift(x)+-self
 def __mul__(self,x):
  x=Jet.lift(x)
  return Jet(self.v*x.v,self.i*x.v+self.v*x.i,self.j*x.v+self.v*x.j,
             self.ij*x.v+self.i*x.j+self.j*x.i+self.v*x.ij)
 __rmul__=__mul__
 def __truediv__(self,x):return self*(1/arb(x))
 def __pow__(self,n):
  out=Jet(1)
  for _ in range(n):out=out*self
  return out
 def diff(self):return Jet(*(x.diff() for x in [self.v,self.i,self.j,self.ij]))

def ej(i,j,k):return Jet(1,k*i,k*j,k*k*i*j)

def local_kernels(p,h,N):
 C=p['C'];B2=p['B2'];lam=p['g'];x=lam*B2;d=p['eps']*C;q=arb(p['q'])
 Fl=1-384*arb.pi()*arb('1e22')*lam;alpha=8*arb.pi()*arb('1e22')/(C*B2)
 Q=mon(0,2)*(p['beta']*(1+2*mon(1)+3*mon(2))+p['u']*mon(2)+p['eps']*mon(4))-p['tau']*mon(3)
 D=lambda f:(Q*f.diff()).diff();R=-Q.diff().diff();V=-mon(4)
 # Complete trace polynomial, with the exact vacuum/product identities preserved.
 du=d/x;eta=lam/(60*arb.pi());dl=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*lam
 fa=4+2/x+p['tau']*C/x;fb=3+1/x+du
 Fn=mon(-2)+1/x-fa*mon(1)+fb*mon(2)-2*du*mon(3)+du*mon(4)
 Lap=lambda f:mon(4)*(Fn*f.diff()).diff()
 rr=4*mon(5)-6*mon(6);rn=-12+du*rr
 tn=(3+1/x+du)*mon(4)-du*(5*mon(5)-arb(9)/2*mon(6))
 cn=fa*mon(3)-2*fb*mon(4)+arb(20)/3*du*mon(5)-5*du*mon(6)
 pp=7*cn**2+arb(22)/3*tn**2-arb(11)/72*(rn**2-144)
 ps=2*(6+du)**2+arb(22)/3*(6+du)/x+arb(22)/3/x**2+22
 S=q*q/(60*arb.pi()*lam*B2*B2)
 tau=S*(-10*mon(4)-20*mon(5)+30*mon(6))+eta*(pp+ps*rr/2)-dl*du*Lap(rr)
 T0=lam*tau
 n=atom('n');u=atom('u');Wdim=d*rr;Rdim=-12*x+Wdim
 D0=-Q.diff().diff()+6*Q.diff()*mon(-1)-12*Q*mon(-2)
 Dv=-2*D(V)-6*D(N)+12*Q*N.diff()*mon(-1)
 Wv=C*mon(2)*(Dv-2*V*D0)-2*N*Rdim
 Wh=-B2*T0-mon(1)*n*Rdim.diff()
 Whv=-8*mon(2)*n*V-2*C*mon(2)*(n+u)*Dv-2*(N+V)*Wh-12*C*mon(2)*Q*n.diff()*N.diff()
 vs=4*n+2*u;vt=4*N+2*V
 sq=alpha*mon(-4)*(2*Wh*Wv+2*Wdim*Whv+2*Wdim*vs*Wv+2*Wdim*vt*Wh+Wdim**2*vs*vt)
 gs=alpha*(8*d*d-16*C*d);sq-=2*gs*u*V
 sqk=sq.strong()
 # Full nonmagnetic Weyl Hessian, before any field derivative is removed.
 phi=Jet(-mon(0,l=1),n,N);bb=Jet(0,u,V)
 LL=Jet(R)-2*D(bb);MM=D(phi);VV=Q*phi.diff()**2
 em=ej(u,V,-2);ep=ej(u,V,2)
 ww=(C*C*em*(phi*LL**2/60+LL*MM/30-MM**2/10+LL*VV/30-7*MM*VV/90-7*VV**2/180)
     -C*phi*LL/18+phi*ep/15-C*VV/18)/(2*arb.pi()*C*B2)
 wk=ww.ij.strong();assert not wk['known'].c
 kM=q*q/(12*arb.pi()*C*B2)
 dk=-(1-Fl)*kM-Fl*(x/C)*eta*ps/2
 kv=11*B2*lam*lam/(120*arb.pi()*C)
 def potential(g,gp,gpp,coef):return {'n':coef*(gpp*N+2*gp*V),'u':coef*(2*gp*N+4*g*V),'known':FP()}
 pk=potential(4*mon(1)-3*mon(2)-1,-4*mon(1)+6*mon(2),4*mon(1)-12*mon(2),dk)
 vk=potential(mon(-4)-1,4*mon(-4),16*mon(-4),kv)
 k=arb(2000);b02=k**10*((1+12*arb.pi()*k**6).sqrt()-1)/6
 Lscale=(B2/b02).log()/2;Feff=4*arb.const_euler()-2*((q+1).digamma()+arb.const_euler())+4*Lscale
 O4=D(D(V))+R*D(V)+D(R*V)+R*R*V/2+2*V
 qk={'n':FP(),'u':q/(6*arb.pi()*B2)*(-D(V)+2*V)+C*Feff/(60*arb.pi()*B2)*O4,'known':FP()}
 parts={'curvature_square':sqk,'nonmagnetic_Weyl':wk,'reference_potential':pk,'vacuum_volume':vk,'LLL_a2_product':qk}
 allk={key:sum((row[key] for row in parts.values()),FP()) for key in ['n','u','known']}
 NUMERICS.append(dict(kind='complete_local_action_coefficients',F_lambda=encoded(Fl),k_M=encoded(kM),delta_k=encoded(dk),kvac=encoded(kv),
  g_sq_product=encoded(gs),F_eff=encoded(Feff),fixed_scale_log=encoded(Lscale),
  combined_kernels={k:v.serial() for k,v in allk.items()},partial_term_counts={n:{k:len(v.c) for k,v in r.items()} for n,r in parts.items()}))
 return allk,Q,O4,Fl,kM,parts

def algebra():
 z,C,x,d,alpha=sp.symbols('z C x d alpha',positive=True)
 n,u,N,V,D0,Dv,Dh,Rdim,T0,B2,Q,np,Np=sp.symbols('n u N V D0 Dv Dh Rdim T0 B2 Q np Np')
 Wh=-B2*T0-z*n*sp.Symbol('Rdim_prime')
 subDh=Wh/(C*z*z)+2*u*D0+2*n*(2/C+D0)
 raw=z*z*(8*n*N+C*(4*(n+u)*(N+V)*D0-2*(n+u)*Dv-2*(N+V)*Dh-12*Q*np*Np))
 target=-8*z*z*n*V-2*C*z*z*(n+u)*Dv-2*(N+V)*Wh-12*C*z*z*Q*np*Np
 check('curvature_mixed_constraint_identity',sp.expand(raw.subs(Dh,subDh)-target)==0)
 b=sp.symbols('b');wp=2-2*C*sp.exp(-2*b)+12*x
 gp=sp.diff(alpha*sp.exp(2*b)*wp**2,b).subs(b,0).subs(C,1+6*x+d).expand()
 check('curvature_square_product_gradient',sp.expand(gp-alpha*(8*d*d-16*(1+6*x+d)*d))==0)
 a,w,Fb,Fn=sp.Function('a')(z),sp.Function('w')(z),sp.Function('Fb')(z),sp.Function('Fn')(z)
 bh=a/2-z*w-z*z*sp.diff(w,z)/2
 check('whole_local_constraint_transpose',sp.expand(z*w*Fn+bh*Fb-a*Fb/2-w*(z*Fn+z*z*sp.diff(Fb,z)/2)+sp.diff(z*z*Fb*w/2,z))==0)

def load_trial(h):
 global POWER
 old=json.loads((root/'evidence/current-schur-scalar-projection.json').read_bytes())['runs'][-1]
 row=next(r for r in old['numerics'] if r['kind']=='actual_global_scalar_inverse')
 POWER=decode(row['fractional_power'])
 return FP({(r['z_power'],r['y_integer_power'],r['fractional_multiplicity'],0):decode(r['coefficient']) for r in row['trial_coefficients']})+mon(4)/2

def correlated_parameters(p,h):
 x=p['g']*p['B2'];old_d=p['eps']*p['C'];dn=old_d
 S=arb(p['q'])**2/(60*arb.pi()*p['g']*p['B2']**2);steps=[]
 for i in range(3):
  du=dn/x;ps=2*(6+du)**2+arb(22)/3*(6+du)/x+arb(22)/3/x**2+22
  nd=x*(-5*S+p['g']*ps/(120*arb.pi()))
  check('product_trace_image_in_original_enclosure_'+str(i),old_d.contains(nd))
  steps.append(encoded(nd));dn=nd
 C=1+6*x+dn;u=1/C;ep=dn/C;beta=(1-u-ep)/6
 v=1+h;tau=v*h*h*(beta*(1+2/v+3/v**2)+u/v**2+ep/v**4)
 out={**p,'C':C,'u':u,'eps':ep,'beta':beta,'gamma':p['kappa']*p['B2'].sqrt()/C,'tau':tau}
 check('same_product_radius_relation',(C-1-6*x-dn).contains(0))
 check('same_product_trace_d_contracts',dn.rad()<old_d.rad()*arb('1e-18'))
 NUMERICS.append(dict(kind='retained_exact_product_parameter_correlations',old_d=encoded(old_d),steps=steps,
  contracted_d=encoded(dn),old_parameters={k:encoded(p[k]) for k in ['C','B2','g','beta','u','eps','gamma','tau']},
  correlated_parameters={k:encoded(out[k]) for k in ['C','B2','g','beta','u','eps','gamma','tau']},
  existing_product_solution_unchanged=True,new_existence_claim=False,normalization_changed=False))
 return out

def J(k,z):
 return ((1-z**(k+1))*(1-z).log()+sum((z**i/i for i in range(1,k+2)),arb(0)))/(k+1)-z**(k+2)/(k+2)

def derivatives(f,n=7):
 out=[f]
 for _ in range(n):out.append(out[-1].diff())
 return out

def joint_range(ds,lo,hi):
 mid=(lo+hi)/2;rad=(hi-lo)/2;c=ds[0].value(mid)
 err=sum((abs(ds[j].value(mid))*rad**j/factorial(j) for j in range(1,7)),arb(0))
 err+=ds[7].absrange(lo,hi)*rad**7/factorial(7)
 return c+arb(0,err.upper())

def outer_integral(poly,zo,weight):
 A=arb('3.768e-12');aloc=arb('4e-28');f=arb('4.34e-11')
 def im(k,l):
  assert k>-1 and l in [0,1]
  out=zo**(k+1)/(k+1)
  return out if l==0 else out*(-zo.log()+arb(1)/(k+1))
 val=arb(0)
 for (k,j,t,l),c in poly.c.items():
  assert k>=0 and l<=1,(k,j,t,l)
  e=j+t*POWER;yb=arb(max(arb(1),((1-zo)**e).upper()))
  if weight=='a':w=A/(2*(1-zo))*im(k+2,l)+aloc*im(k,l)
  else:w=weight*im(k+1,l)+A/(1-zo)*im(k+2,l)+2*(f+aloc)*im(k,l)
  val+=abs(c)*yb*w
 return val

def cold_moment(e,y0,yc,log_weight=False):
 a=e+1
 if a.is_zero():
  i=(yc/y0).log()
  return i if not log_weight else (1+arb('.01').log())*i-(yc.log()**2-y0.log()**2)/2
 if not log_weight:return (yc**a-y0**a)/a
 F=lambda y:y**a/a*(1+(arb('.01')/y).log()+1/a)
 return F(yc)-F(y0)

def cold_bound(poly,y0,yc,with_log):
 zc=1-yc;out=arb(0)
 for (k,j,t,l),c in poly.c.items():
  assert l<=1
  zp=arb(max(arb(1),(zc**k).upper()))
  mom=cold_moment(j+t*POWER+1+l,y0,yc,with_log)
  assert mom>=0,(k,j,t,l,str(mom))
  out+=abs(c)*zp/zc**(1+l)*mom
 return out

def integrate_local(kernels,h,zh,Fl):
 A=arb('3.768e-12');aloc=arb('4e-28');f=arb('4.34e-11');s0=arb('3.33e-10')*h
 zo=arb(2)**-40;yc=arb(2)**-49;zc=1-yc;y0=h/(1+h)
 ka=kernels['u']/2;kw=kernels['n']+mon(1)*kernels['u'].diff()/2
 known=kernels['known'].integral(y0)
 check('known_curvature_integral_roundoff_budget',2*known.rad()/Fl<arb('1e-44'),known/Fl)
 da=derivatives(ka);dw=derivatives(kw)
 outer=outer_integral(ka,zo,'a')+outer_integral(kw,zo,s0)
 def cell(lo,hi):
  ra=abs(joint_range(da,lo,hi));rw=abs(joint_range(dw,lo,hi))
  wa=A*(J(0,hi)-J(0,lo))+aloc*(hi-lo)
  ww=s0*(hi*hi-lo*lo)/2+2*A*(J(0,hi)-J(1,hi)-J(0,lo)+J(1,lo))+2*(f+aloc)*(hi-lo-(hi*hi-lo*lo)/2)
  assert wa>0 and ww>0
  lower=arb(max(arb(0),ra.lower()))*wa+arb(max(arb(0),rw.lower()))*ww
  upper=arb(ra.upper())*wa+arb(rw.upper())*ww
  return dict(lo=lo,hi=hi,lower=lower,upper=upper,priority=float((upper-lower).upper()))
 heap=[];serial=0;lower=arb(0);upper=arb(outer.upper())
 def add(row):
  nonlocal serial,lower,upper
  serial+=1;heapq.heappush(heap,(-row['priority'],serial,row));lower+=row['lower'];upper+=row['upper']
 for k in range(40,1,-1):add(cell(arb(2)**-k,arb(2)**(-k+1)))
 for k in range(2,50):add(cell(1-arb(2)**(-k+1),1-arb(2)**-k))
 while (upper-lower).upper()>arb('5e-44'):
  assert len(heap)<30000,'local complete-cell refinement limit'
  _,_,row=heapq.heappop(heap);lower-=row['lower'];upper-=row['upper'];mid=(row['lo']+row['hi'])/2
  add(cell(row['lo'],mid));add(cell(mid,row['hi']))
  if len(heap)%2000==0:print(json.dumps(dict(event='local_contact_refinement',cells=len(heap),enclosure_width=str(upper-lower),upper=str(upper))),flush=True)
 rows=sorted((r for _,_,r in heap),key=lambda r:float(r['lo']))
 check('local_complete_bulk_cover',rows[0]['lo']==zo and rows[-1]['hi']==zc and all(a['hi']==b['lo'] for a,b in zip(rows,rows[1:])))
 check('local_positive_integral_enclosure_precision',upper-lower<arb('5e-44'))
 ell=-yc.log()-zc;wb=s0+2*yc/zc*(A*ell+f+aloc)
 endpoint=abs(zc*zc*kernels['u'].value(zc)/2)*wb
 cold=arb('1e-6')*cold_bound(kernels['u'],y0,yc,True)+arb('1.2e-8')*cold_bound(kernels['n'],y0,yc,False)
 allowance=(abs(known)+upper+endpoint+cold)/Fl
 NUMERICS.append(dict(kind='whole_domain_local_projection_enclosure',
  known_signed_integral=encoded(known/Fl),unknown_actual_pairing_allowance=encoded((upper+endpoint+cold)/Fl),
  total_local_absolute_allowance=encoded(allowance),local_preassigned_budget='1e-41',meets_local_budget=bool(allowance<arb('1e-41')),
  bulk_absolute_integral=encoded(hull(lower,upper)),analytic_outer_allowance=encoded(outer),retained_cut_endpoint=encoded(endpoint),
  complete_cold_allowance=encoded(cold),physical_y_bolt=encoded(y0),outer_cut=encoded(zo),cold_cut=encoded(zc),
  cells=[dict(lo=encoded(r['lo']),hi=encoded(r['hi']),lower=encoded(r['lower']),upper=encoded(r['upper'])) for r in rows],
  full_actual_h1_signed_value_evaluated=False))
 print(json.dumps(dict(event='whole_local_contact_enclosure',cells=len(rows),known=str(known/Fl),unknown_pairing_upper=str((upper+endpoint+cold)/Fl),
  total_upper=str(allowance),cold=str(cold),endpoint=str(endpoint),local_budget_passed=bool(allowance<arb('1e-41')))),flush=True)
 return allowance

def absolute_poly_integral(f):
 out=arb(0)
 for (k,j,t,l),c in f.c.items():
  assert k>=0 and j>=0 and t==0 and l==0,(k,j,t,l)
  out+=abs(c)*arb(factorial(k)*factorial(j))/factorial(k+j+1)
 return out

def angular_tail(p,Q,O4,Fl):
 B2=p['B2'];C=p['C'];q=arb(p['q']);H=arb('1.2e-8');t0=arb(1)/1000
 old=json.loads((root/'evidence/current-disk-heat-propagation.json').read_bytes())['runs'][-1]
 row=next(r for r in old['numerics'] if r['kind']=='complete_actual_domain_propagation_budget')
 point_error=decode(row['first_actual']['base'])+decode(row['higher_Volterra'])
 check('positive_time_trace_anchor_error_guard',point_error<100)
 ta=t0/4;trace_guard=2*(arb('.19')/ta+100*ta**3)
 check('actual_full_thermal_trace_anchor',trace_guard<1521)
 V=-mon(4);Dv=(Q*V.diff()).diff()
 S1=H*absolute_poly_integral(Dv-2*V)/(12*arb.pi())
 S2=H*absolute_poly_integral(O4)/(120*arb.pi())
 I=arb('6.704e-12');B0=arb('2.3454197e-13')*t0**2/4+arb('1.3212e-5')*t0**3/6
 Ch=2*arb(2).sqrt()*arb('2.25')*1521*H/arb(1).exp();out=[]
 check('complete_angular_tail_monotonicity',2*C<t0*q*q)
 for N in [32,64,80,128]:
  lam=arb(N+1)*(N+1+q)/C
  theta=(-lam*t0).exp()*(2*(q+2*N+2)+2*C/t0)
  err=2*theta/(B2*Fl)*(B0+Ch*t0.exp()*lam/(lam-1)+I/(arb.pi()*t0)+I/(6*arb.pi())+S1+S2*(t0+1/lam))
  out.append(dict(N=N,first_omitted_lambda=encoded(lam),complete_theta_tail=encoded(theta),complete_positive_time_tail=encoded(err),
    fits_existing_reserve=bool(err<arb('2.5e-43')),fits_tenth_of_reserve=bool(err<arb('2.5e-44'))))
 check('complete_positive_time_tail_at_80',decode(out[2]['complete_positive_time_tail'])<arb('2.5e-44'))
 NUMERICS.append(dict(kind='actual_complete_positive_time_angular_tail',trace_anchor_time=encoded(ta),
  actual_trace_anchor_upper=encoded(trace_guard),Duhamel_trace_coefficient=encoded(Ch),s1_upper=encoded(S1),s2_upper=encoded(S2),
  finite_time_boundary_upper=encoded(B0),rows=out,low_angular_signed_projection_evaluated=False,
  all_physical_thermal_modes_retained=True,full_angular_tail_included=True))
 print(json.dumps(dict(event='complete_positive_time_angular_tail',upper_by_N={r['N']:str(decode(r['complete_positive_time_tail'])) for r in out})),flush=True)

def execute():
 p,h,zh,env=gy.actual_parameters();NUMERICS.extend(gy.NUMERICS)
 p=correlated_parameters(p,h)
 N=load_trial(h);kernels,Q,O4,Fl,kM,parts=local_kernels(p,h,N)
 yc=arb(2)**-49;zc=1-yc
 # Global guards avoid any change of the nonlocal operator at zc.
 kmcold=kM/Fl*arb('1.2e-8')*yc*(16*arb('1.01')+24+24*arb('1.01')+36)
 check('whole_domain_magnetic_cold_fits_unchanged_reserve',kmcold<arb('2e-43'),kmcold)
 A=arb('3.768e-12');a=arb('4e-28');f=arb('4.34e-11');sr=arb('3.33e-10')*h
 wb=sr+2*yc/zc*(A*(-yc.log()-zc)+f+a)
 pair=A*(arb(29)/200+4*arb(103)/4200)+4*f/30+a*(arb(1)/10+arb(4)/30)+sr/3+zc**6*wb/2+arb('1.2e-8')*yc
 check('global_product_pairing_includes_cold_tip',pair<arb('6.704e-12'),pair)
 NUMERICS.append(dict(kind='whole_domain_projection_split',magnetic_cold_upper=encoded(kmcold),
  full_pairing_upper=encoded(pair),rest_target='1.287e-40',old_complete_cold_not_added_again=True,
  discontinuous_metric_cutoff_used=False,normalization_changed=False))
 angular_tail(p,Q,O4,Fl)
 integrate_local(kernels,h,zh,Fl)

def main():
 ctx.prec=512
 old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 review=json.loads((root/OWN[2]).read_bytes())
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
  {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=bh.prior.scalar.prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  dependencies=deps,invocation=sys.argv,frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
 def save():out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 save()
 try:
  check('previous_scoped_execution_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('interval_library_pinned',bh.prior.scalar.prev.prior.flint.__version__=='0.9.0')
  algebra();execute()
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=bh.prior.scalar.prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  continuum_argument_scope='Ordinary complete-action, actual-domain smearing and Schatten/Duhamel argument; checks only support the displayed identities and inequalities.',
  full_local_contact_enclosure_executed=status=='passed',complete_signed_projection_evaluated=False,
  actual_controlled_solution=False,nonlinear_action_domain_proved=False,goal_complete=False)
 save()
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
