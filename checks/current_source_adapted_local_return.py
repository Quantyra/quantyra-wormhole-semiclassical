"""Whole-domain two-row return enclosure, with the fractional flow entry separate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import argparse, json, sys, traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
from flint import acb
from current_source_adapted_geodesic_first_attempt import result,constants
from current_localized_first_residual_range_test import pa,ps,pm,pv,pd,mag

ROOT=Path(__file__).resolve().parents[1]
STEM='current-source-adapted-local-return'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_source_adapted_local_return.py',
 'checks/current_quartic_ritz_rank_test.py','checks/current_source_adapted_geodesic_first_attempt.py',
 'checks/current_localized_first_residual_range_test.py','checks/current_localized_short_time_trial.py',
 'notes/current-finite-local-action-assembly.md','notes/current-canonical-heat-restoration-assembly.md',
 'notes/current-geodesic-two-sided-remainder.md','notes/current-source-adapted-geodesic-first-attempt.md',
 'evidence/current-source-adapted-geodesic-first-attempt.json',
 'evidence/current-source-adapted-geodesic-first-attempt-verification.json',
 'evidence/current-global-green-trial-enclosure.json','evidence/current-global-mass-conversion.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p):return sha256(p.read_bytes()).hexdigest()
def up(x):return arb(abs(x).upper())
def minimum(*x):return arb(min(up(v).upper() for v in x))
def maximum(*x):return arb(max(up(v).upper() for v in x))

class J:
 """Taylor coefficients, with explicit loss of known order on differentiation."""
 def __init__(self,c):
  assert not any(isinstance(v,float) for v in c),'Binary floats are not certified coefficients'
  self.c=[arb(v) if isinstance(v,int) else v for v in c]
 @staticmethod
 def cast(x,n):return x if isinstance(x,J) else J([x]+[0]*(n-1))
 def __add__(self,o):
  if isinstance(o,V):return NotImplemented
  o=J.cast(o,len(self.c));return J([a+b for a,b in zip(self.c,o.c)])
 __radd__=__add__
 def __neg__(self):return J([-a for a in self.c])
 def __sub__(self,o):return self+-J.cast(o,len(self.c))
 def __rsub__(self,o):return -self+o
 def __mul__(self,o):
  if isinstance(o,V):return NotImplemented
  o=J.cast(o,len(self.c));n=min(len(self.c),len(o.c))
  return J([sum((self.c[j]*o.c[i-j] for j in range(i+1)),0) for i in range(n)])
 __rmul__=__mul__
 def inv(self):
  r=[1/self.c[0]]
  for i in range(1,len(self.c)):r.append(-sum((self.c[j]*r[i-j] for j in range(1,i+1)),0)/self.c[0])
  return J(r)
 def __truediv__(self,o):return self*J.cast(o,len(self.c)).inv()
 def __rtruediv__(self,o):return self.inv()*o
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  r=J.cast(1,len(self.c))
  for _ in range(n):r=r*self
  return r
 def d(self):
  assert len(self.c)>1,'Insufficient coefficient differentiation order'
  return J([i*self.c[i] for i in range(1,len(self.c))])
 def log(self):
  derivative=self.d()/self
  return J([self.c[0].log()]+[derivative.c[i-1]/i for i in range(1,len(self.c))])

NB=8
class V:
 """Exact linear expression in r,hr,f,hf,g,hg,k,W."""
 def __init__(self,c):self.c=list(c)
 def __add__(self,o):
  if o==0:return self
  assert isinstance(o,V);return V([a+b for a,b in zip(self.c,o.c)])
 __radd__=__add__
 def __neg__(self):return V([-a for a in self.c])
 def __sub__(self,o):return self+-o
 def __mul__(self,o):return V([a*o for a in self.c])
 __rmul__=__mul__
 def __truediv__(self,o):return self*(arb(1)/o)

def rows(p,c,bs,x,phi_value=None):
 """Signed strong rows on a whole real cell or a complex coefficient arc."""
 L,a,ga,beta,C,B2=[p[k] for k in ['L','a','gamma','beta','C','B2']]
 z=L-ga*x;q2=2*x+a*x*x
 # Preserve the exact quadratic cancellation before interval evaluation.
 rq=1+sum((bs[j]*ga**(j-2)*x**(j-1)/(2+a*x) for j in range(3,7)),0)
 q=q2*rq;qx=q.d();qxx=qx.d();D=z*qxx.d()+3*ga*qxx
 def unit(i):return V([J.cast(int(j==i),len(x.c)) for j in range(NB)])
 r,hr,f,hf,g,hg,k,W=[unit(i) for i in range(NB)]
 hb=[hr,q2*(c['muY']*r-f),hf,-q2*g,hg,-q2*k,None,q2*D*r]
 db=[c['muY']*r-f,c['muY']*hr-hf+q2.d()*(c['muY']*r-f),
     -g,-hg-q2.d()*g,-k,None,None,D*hr+(q2*D.d()+q2.d()*D)*r]
 def H(v):
  out=V([J.cast(0,len(x.c)) for _ in range(NB)])
  for i,A in enumerate(v.c):
   if all(t==0 for t in A.c):continue
   assert hb[i] is not None,('unsupported first derivative',i)
   out=out+unit(i)*(q2*A.d())+hb[i]*A
  return out
 def Delta(v):
  out=V([J.cast(0,len(x.c)) for _ in range(NB)])
  for i,A in enumerate(v.c):
   if all(t==0 for t in A.c):continue
   assert db[i] is not None and hb[i] is not None,('unsupported second derivative',i)
   out=out+unit(i)*(q2*A.d().d()+q2.d()*A.d())+hb[i]*(2*A.d())+db[i]*A
  return out
 def Da(v):return rq*Delta(v)+rq.d()*H(v)
 M=z*qxx+4*ga*qx
 nu=z**4/(3*L**5)*(-z*rq*f+(z*rq*c['muY']-M)*r-(7*ga*rq-z*rq.d())*hr+W)
 u=z**4/L**5*(z*rq*f+(2*M/3-z*rq*c['muY'])*r+(5*ga*rq-z*rq.d())*hr-2*W/3)
 dn,du,hn=Da(nu),Da(u),H(nu)
 L0=-qxx;M0=ga*qx/z+ga**2*q/z**2;N0=ga**2*q/z**2
 D0=L0-6*M0-6*N0;lv=-2*du;mv=dn;nv=2*ga*rq/z*hn
 dv=lv-6*mv-6*nv
 l0=arb(2000)**-10*B2/C;alpha=8*arb.pi()*arb('1e22')/(C*B2)
 zmag=arb.pi()*arb(2000)**16/(C*B2)
 fnE=(2*D0*nu+dv)/z**2+4*(nu+u)/(C*z**2)+12*l0*(4*nu+2*u)/z**4
 fuE=-2*Da(nu/z**2)+4*(nu+u)/(C*z**2)+12*l0*(2*nu+u)/z**4-4*zmag*u
 w=C*(dv-2*D0*u);U0=2+C*D0
 def divphi(v):return ga*(rq/z*H(v)+(q/z).d()*v)
 def divnu(A):return A*rq*Delta(nu)+(A*rq).d()*hn
 fnR=-12*alpha*C*(Da(w)-2*divphi(w)-2*divnu(U0))
 fuR=alpha*(16*u+4*C*C*D0**2*u-4*C*C*D0*dv-4*C*Da(w))
 modeln=c['fn_r']*r+c['fn_f']*f+c['fn_g']*g
 modelu=c['fu_r']*r+c['fu_f']*f+c['fu_g']*g
 classicaln=-(fnE+fnR)/beta-modeln
 classicalu=-(fuE+fuR)/beta-modelu
 phi=-z.log()
 if phi_value is not None:phi.c[0]=phi_value
 km=arb(1)/15-arb(2000)**8/3
 A0=L0/30-7*M0/90-7*N0/90
 av=lv/30-7*mv/90-7*nv/90-2*A0*u
 fnW=C*C*((2*L0*lv-2*L0**2*u)/60+Da(lv-2*L0*u)/30
       +Da(-mv/5-7*nv/90-2*(-M0/5-7*N0/90)*u)-2*divphi(av)-2*divnu(A0))
 fnW=fnW-C*lv/18+2*km*u+C*mv/9
 B0=phi*L0**2/60+L0*M0/30-M0**2/10+L0*N0/30-7*M0*N0/90-7*N0**2/180
 bv=L0**2*nu/60+phi*L0*lv/30+(M0*lv+L0*mv)/30-M0*mv/5
 bv=bv+(N0*lv+L0*nv)/30-7*(N0*mv+M0*nv)/90-7*N0*nv/90
 fuW=C*C*(4*B0*u-2*bv-2*Da((L0*nu+phi*lv+mv+nv-2*(phi*L0+M0+N0)*u)/30))
 fuW=fuW+C*dn/9+2*km*(nu+2*phi*u)
 fnW=-fnW/(2*arb.pi()*C*B2*beta);fuW=-fuW/(2*arb.pi()*C*B2*beta)
 # Retain all common local restorations in their known coefficient intervals.
 c0=2*(604*arb(2000)**8+(1+2*arb('1e-8'))*(arb(2000)**8/6-arb(1)/30))/(arb.pi()*C*B2)
 b2=arb(0,100)*C/(480*arb.pi()*B2)
 P4=4*L0**2*u+8*L0*du+8*Da(L0*u)+8*Da(du)
 fuW=fuW-arb(0,up(c0))/beta*u-b2/beta*P4+arb(2000)**4/(6*arb.pi()*B2*beta)*du
 return [[v.c[i].c[0]*z.c[0]**2 for i in range(NB)] for v in [classicaln,classicalu,fnW,fuW]]

def model_rows(p,c):
 a,L,C,B2,beta=[p[k] for k in ['a','L','C','B2','beta']]
 l0=arb(2000)**-10*B2/C;Z=arb.pi()*arb(2000)**16/(C*B2)
 al=8*arb.pi()*arb('1e22')/(C*B2)
 nn=[-4*a/L**2+4/(C*L**2)+48*l0/L**4,6/L**2+al*(-48*C+48*a*C*C),72*al*C*C]
 nu=[4/(C*L**2)+24*l0/L**4,2/L**2+48*al*a*C*C,24*al*C*C]
 uu=[4/(C*L**2)+12*l0/L**4-4*Z+al*(16+16*a*a*C*C),32*al*a*C*C,8*al*C*C]
 def compose(n,m):
  assert (m[2]-n[2]/3).contains(0)
  return [-(m[0]-n[0]/3-(2*a/9)*(n[1]-c['muY']*n[2]))/beta,
          -(m[1]-n[1]/3-(2*a/9)*n[2])/beta,
          (2*a/9)*(n[0]-c['muY']*n[1]+c['muY']**2*n[2])/beta]
 for name,n,m in [('fn',nn,nu),('fu',nu,uu)]:
  c[name+'_f'],c[name+'_g'],c[name+'_r']=compose(n,m)

def profile_terms(p,c,name):
 a=p['a'];ll=a*p['ell']*(p['ell']+1);ln=c['large_mass_lambda']
 poles=[c['muT'],c['muS'],ll];num=[c['K0']]
 if name in ['r','rf','fn','fu']:poles.append(c['muY'])
 if name in ['fn','fu','rf']:
  label='fu' if name=='rf' else name
  num=ps(pa([c[label+'_r']],pm([c[label+'_f'],c[label+'_g']],[c['muY'],arb(1)])),c['K0'])
 if name=='rf':poles.append(2*a)
 if name=='g':num=pm(num,[arb(0),arb(1)])
 if name=='k':num=pm(num,[arb(0),arb(0),arb(1)])
 return [(num,poles),(ps(num,(c['Cv']-1)*ln),poles+[ln])]

def profile_components(p,c,name):
 """Grouped residue contour avoids division by almost coincident masses."""
 a=p['a'];terms=profile_terms(p,c,name);rho=arb('.25');center=arb(2)
 weight=arb(0);isopoles=[]
 for num,poles in terms:
  numerator=sum((up(b)*arb('2.25')**i for i,b in enumerate(num)),arb(0));den=arb(1)
  for mu in poles:
   dist=abs(mu-center)
   if dist<arb('.01'):den*=rho-dist
   else:
    assert dist>rho;den*=dist-rho
    if not any(mu.overlaps(v) for v in isopoles):isopoles.append(mu)
  assert den>0;weight+=rho*numerator/den
 components=[(arb('.9'),up(weight),arb(2))]
 for mu in isopoles:
  residue=arb(0)
  for num,poles in terms:
   matches=[i for i,v in enumerate(poles) if mu.overlaps(v)]
   if not matches:continue
   assert len(matches)==1
   den=arb(1)
   for j,v in enumerate(poles):
    if j!=matches[0]:den*=v-mu
   assert not den.contains(0);residue+=pv(num,-mu)/den
  k=((1+4*mu/a).sqrt()-1)/2
  components.append((k,up(residue),arb(1)))
 return components

def profile_bounds(p,c,comps,s0,s1,kind='norm'):
 a,ga,At=p['a'],p['gamma'],p['At'];width=s1-s0
 assert width>0
 X0=a*s0/ga;X1=a*s1/ga;ans=arb(0)
 for k,w,fluxfactor in comps:
  N=(ga/a*acb(k+1).polygamma(1).real/(2*k+1)).sqrt()
  if kind=='norm':
   if X1<=1:
    def primitive(x):
     if x==0:return arb(0)
     v=(3/x).log();return x*(v*v+2*v+2)/4
    bound=(ga/a*(primitive(X1)-primitive(X0))).sqrt()
   elif X0>0:
    bound=width.sqrt()*(1+2/X0).log()/2
    if X0>=1:bound=minimum(bound,width.sqrt()*arb.pi()*X0**(-k-1))
   else:bound=N
   ans+=w*minimum(N,bound)
  elif kind=='flux':
   bound=width.sqrt()
   if X0>=1:bound=minimum(bound,width.sqrt()*arb.pi()*(k+1)*X0**(-k))
   ans+=w*fluxfactor*bound
  elif kind=='tail_integral':
   bound=ga/(a*k*(k+1))
   if X0>=1:bound=minimum(bound,arb.pi()*ga/(a*k)*X0**(-k))
   ans+=w*bound
  else:raise ValueError(kind)
 return up(At*ans)

def moments(p,c,comps):
 a,ga,L,At=[p[k] for k in ['a','gamma','L','At']]
 X=a*L/ga;B=(arb(3).log()**2/3+2*arb(3).log()/9+arb(2)/27)/4
 N=M=H=arb(0)
 for k,w,ff in comps:
  tint=1/(2*k-1) if k>arb('.5') else (X**(1-2*k)-1)/(1-2*k)
  N+=w*(ga/a*acb(k+1).polygamma(1).real/(2*k+1)).sqrt()
  M+=w*((ga/a)**3*(B+arb.pi()**2*tint)).sqrt()
  H+=w*ff*ga*(ga/a*(1+arb.pi()**2*(k+1)**2*tint)).sqrt()
 return [up(At*v) for v in [N,M,H]]

def execute():
 ctx.prec=512
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 p['At']=dec(result('current-global-green-trial-enclosure')['constants']['A_t'])
 mp={k:dec(v) for k,v in result('current-global-mass-conversion')['parameters'].items()}
 c=constants('current-source-adapted-geodesic-first-attempt');model_rows(p,c)
 L,a,ga,beta=[p[k] for k in ['L','a','gamma','beta']]
 u,eps,tau=[mp[k] for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+eps,-2*eps,eps]
 bs=[sum((Q[j]*comb(j,i)*L**(j-i)*(-1)**i for j in range(i,7)),arb(0)) for i in range(7)]
 assert bs[0].contains(0) and bs[1].overlaps(2*ga) and bs[2].overlaps(a)
 bs[0]=arb(0);bs[1]=2*ga;bs[2]=a
 assert a>arb('.999') and a<arb('1.001') and bs[3]>arb('-1.34') and bs[3]<arb('-1.32')
 assert bs[4]>arb('.49') and bs[4]<arb('.51') and abs(bs[5])+abs(bs[6])<arb('1e-25')
 assert ga<arb('1e-30') and L<=1 and L>arb('.99')
 # Lower guards for Pxi and QPxi, after the exact polynomial cancellations.
 Pguard=arb('1.5')*beta-arb('4.5')*abs(u)-20*(abs(tau)+abs(eps))
 QPguard=12*beta-20*(abs(u)+abs(eps))
 assert Pguard>arb('.24') and QPguard>1
 Cxi=arb('1.01');cy=1+2/arb('4.5').sqrt()
 CU=Cxi*(arb('.5')+cy/arb('4.5').sqrt());CV=Cxi*(1+cy/arb('4.5').sqrt()+cy/2)
 # The full polynomial also bounds Q, Q''' and the Volterra coefficient.
 assert beta<arb('.17')
 M3=sum((j*(j-1)*(j-2)*up(Q[j]) for j in range(3,7)),arb(0))
 assert M3<1000
 Dmax=sum((j*(j-1)*abs(5-j)*up(Q[j]) for j in range(2,7)),arb(0))
 names=['r','f','g','k','fn','fu','rf'];profiles={n:profile_components(p,c,n) for n in names}
 mm={n:moments(p,c,profiles[n]) for n in ['fn','fu','rf']}
 def wf(n):
  N,M,H=mm[n];return (2*ga).sqrt()*(N*M).sqrt()+a.sqrt()*M
 # Complete spectral estimate for sqrt(Q2) rf_s, retaining every rational factor.
 m=arb(min(arb(1).lower(),c['muT'].lower(),c['muS'].lower(),(a*p['ell']*(p['ell']+1)).lower()))
 poly=pa([c['fu_r']],pm([c['fu_f'],c['fu_g']],[c['muY'],arb(1)]))
 ratio=sum((up(v)/m**(2-i) for i,v in enumerate(poly)),arb(0))
 Crf=c['K0']*c['Cv']*ratio
 Grf=p['At']*ga.sqrt()*Crf/(arb(40).sqrt()*m*m)
 CN=maximum(arb(8),arb(4000)/3+8*a+2000*L)
 Hrf=mm['rf'][2]
 Eadj=2*Cxi*Hrf+cy*Cxi*(arb(2).sqrt()*(wf('fn')+2*(2*a*wf('rf')+wf('fu')))
     +2*arb('.17').sqrt()*Hrf+CN/(2*arb('.1').sqrt())*(2*ga*Grf+2*Hrf/a.sqrt()))
 scale=c['original_model_fractional_scale']
 def basis_bounds(s0,s1):
  ns={n:minimum(profile_bounds(p,c,profiles[n],s0,s1),scale) for n in ['r','f','g','k']}
  qmax=2*s1/ga+a*(s1/ga)**2
  hs={n:minimum(profile_bounds(p,c,profiles[n],s0,s1,'flux'),qmax.sqrt()*scale) for n in ['r','f','g']}
  W=Dmax*(s1-s0).sqrt()*profile_bounds(p,c,profiles['r'],s0,s1,'tail_integral')
  return [ns['r'],hs['r'],ns['f'],hs['f'],ns['g'],hs['g'],ns['k'],W]
 # Cover the entire near/middle domain by whole interval cells, no sampled norms.
 smax=L-arb('.01');edges=[arb(0),ga]
 while True:
  last=edges[-1]
  nxt=arb(min((2*last).mid(),(last+(L-last)/10).mid()))
  if not nxt<smax:break
  edges.append(nxt)
 edges.append(smax)
 totals=[arb(0)]*4;cells=[]
 for n,(s0,s1) in enumerate(zip(edges,edges[1:])):
  x0=s0/ga;x1=s1/ga;x=J([arb((x0+x1)/2,(x1-x0)/2),arb(1)]+[arb(0)]*8)
  rr=rows(p,c,bs,x);bb=basis_bounds(s0,s1)
  costs=[sum((up(t)*b for t,b in zip(row,bb)),arb(0)) for row in rr]
  assert all(v.is_finite() for v in costs),(n,costs)
  totals=[v+w*w for v,w in zip(totals,costs)]
  cells.append(dict(s0=enc(s0),s1=enc(s1),row_norm_upper=[enc(v) for v in costs]))
  if n%16==0:print(json.dumps({'cell':n,'of':len(edges)-1,'costs':[str(v) for v in costs]}),flush=True)
 # Entire outer coefficient circle, split into 64 interval arcs.
 outerA=[[arb(0)]*NB for _ in range(4)];outerB=[[arb(0)]*NB for _ in range(4)]
 for j in range(64):
  theta=arb((arb(2)*j+1)*arb.pi()/64,arb.pi()/64)
  zz=acb(theta.cos(),theta.sin())*arb('.02')
  xx=J([(L-zz)/ga,acb(1)]+[acb(0)]*8)
  r0=rows(p,c,bs,xx,acb(0));r1=rows(p,c,bs,xx,acb(1))
  for h in range(4):
   for i in range(NB):
    outerA[h][i]=maximum(outerA[h][i],abs(r0[h][i]))
    outerB[h][i]=maximum(outerB[h][i],abs(r1[h][i]-r0[h][i])) if h>=2 else arb(0)
 bb=basis_bounds(L-arb('.01'),L)
 # On this short outer interval all profiles have a uniform monotone tail bound.
 # Use sup values from the left s endpoint, not an L2-log product inequality.
 s0=L-arb('.01');X=a*s0/ga
 point=[]
 for name,flux in [('r',False),('r',True),('f',False),('f',True),('g',False),('g',True),('k',False)]:
  val=arb(0)
  for kk,ww,ff in profiles[name]:
   val+=ww*(ff*arb.pi()*(kk+1)*X**(-kk) if flux else arb.pi()*X**(-kk-1))
  point.append(p['At']*val)
 point.append(Dmax*profile_bounds(p,c,profiles['r'],s0,L,'tail_integral'))
 b=arb('.01');lognorm=(b*(b.log()**2-2*b.log()+2)).sqrt()
 outer=[]
 for h in range(4):
  val=sum(((outerA[h][i]*b.sqrt()+outerB[h][i]*lognorm)*point[i] for i in range(NB)),arb(0))
  assert val.is_finite();outer.append(val);totals[h]+=val*val
 norms=[t.sqrt() for t in totals]
 Eclassical=CU*norms[0]+CV*norms[1];Equantlocal=CU*norms[2]+CV*norms[3]
 Ereturn=Eclassical+Equantlocal+Eadj
 target=arb('4.313e32');remaining=target-Ereturn
 kout=c['actual_reference_fractional_output_coefficient']
 vals=dict(CU=CU,CV=CV,Pxi_positive_guard=Pguard,QPxi_guard=QPguard,
  Dmax=Dmax,complete_classical_two_row_difference=Eclassical,complete_local_quantum_rows=Equantlocal,
  complete_adjoint_York_return=Eadj,complete_local_return=Ereturn,
  local_return_over_remaining_account=Ereturn/target,
  account_after_local_NOT_physical_margin=remaining,
  first_stage_known_subtotal=c['known_first_stage_subtotal']+Ereturn,
  local_original_bolt_estimator_debit=c['first_residual_to_bolt_multiplier']*Ereturn)
 if remaining>0:vals['conditional_actual_flow_norm_ceiling_APPLICATION_ZERO']=remaining/kout
 return ('complete_local_return_enclosed_fractional_entry_unknown' if Ereturn<target else
         'complete_local_return_estimator_excessive'),dict(constants={k:enc(v) for k,v in vals.items()},
  cell_count=len(cells),cells=cells,outer_row_norm_upper=[enc(v) for v in outer],
  model_row_coefficients={k:enc(v) for k,v in c.items() if k.startswith(('fn_','fu_'))},
  profile_components={n:[dict(degree=enc(k),weight=enc(w),flux_factor=enc(ff)) for k,w,ff in v] for n,v in profiles.items()},
  profile_moments={n:[enc(v) for v in row] for n,row in mm.items()},
  actual_fractional_flow_norm=None,actual_application_error=None,complete_first_reference_residual=None,
  scope=dict(full_local_return=True,full_reference_inverse=False,later_actual_quantum_applications=False,
   original_admission=False,nonlinear_existence=False,compatibility=False,physical_matching=False,
   validity_noise=False,all_sector_stability=False,construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 path=ROOT/'evidence'/(args.receipt+'.json');assert not path.exists(),'Preserve immutable receipts'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());assert digest(Path(sys.executable))==rt['python_sha256']
 pkg=rt['packages']['flint'];assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['status'],rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 path.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if rec['status']=='failed':raise SystemExit(1)
