"""Actual complete localized short-time column, with full analytic tails."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import json,sys,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
import current_disk_heat_propagation as d
import current_localized_short_time_trial as l
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-quantum-source-test'
INPUTS=['notes/'+STEM+'.md','checks/current_localized_quantum_source_test.py',
 'checks/current_quartic_ritz_rank_test.py','checks/current_disk_heat_propagation.py',
 'checks/current_localized_short_time_trial.py','notes/current-disk-heat-propagation.md',
 'notes/current-local-boundary-heat-remainder.md','notes/current-localized-kernel-return-attempt.md',
 'notes/current-localized-boundary-layer-test.md','evidence/current-disk-heat-propagation.json',
 'evidence/current-global-green-trial-enclosure.json','evidence/current-global-mass-conversion.json',
 'evidence/current-localized-short-time-trial.json','evidence/current-positive-time-quantum-operator.json',
 'evidence/current-localized-boundary-layer-test.json','evidence/current-localized-kernel-return-attempt.json',
 'evidence/current-first-quantum-solve-attempt.json','evidence/current-first-inverse-source-screen.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
CELLS=[]
def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(s):
 r=json.loads((ROOT/'evidence'/(s+'.json')).read_bytes());r=r.get('runs',[r])[-1]
 assert r['status']!='failed';return r['result']
def aa(x):return arb(abs(x).upper())
def up(x):return arb(x.upper())
def mx(*xs):return arb(max(x.upper() for x in xs))
def read(xs):return {(x['power'],x['log_power']):dec(x['coefficient']) for x in xs}
def E(p):return l.mul(l.mono(1),l.diff(p))
def sup(p,wlo,whi,shift=0):
 lg=aa(wlo.log())
 return sum((aa(c)*(wlo if n-shift<0 else whi)**(n-shift)*lg**k for (n,k),c in p.items()),arb(0))
def peval(cs,x):return d.pv(cs,x)
def pair(x):return [dec(x['inf']),dec(x['row'])]
def cv(p,f):return p[0]*f[1]/32+7*p[1]*f[0]/24

def Udata(v,wlo,whi):
 terms=[];cs=[arb(0)]*(d.N+3)
 for (n,k),coef in v.items():
  assert k<=1
  a=aa(coef)*(wlo if n-3<0 else whi)**(n-3);lg=aa(wlo.log())
  ts=[]
  for j in range(d.N+3):
   ts.append(arb(comb(n+j-1,j)) if n>0 else arb(int(j==0)) if n==0 else arb(comb(-n,j) if j<=-n else 0))
  for j in range(d.N+3):
   z=ts[j] if not k else lg*ts[j]+sum((ts[j-i]/i for i in range(1,j+1)),arb(0))
   cs[j]+=a*z
  terms.append((n,k,a,lg))
 def norms(x):
  out=[arb(0)]*3
  for n,k,a,lg in terms:
   if n>=0:
    t=(1-x)**(-n);t1=n*t/(1-x);t2=n*(n+1)*t/(1-x)**2
   else:
    m=-n;t=(1+x)**m;t1=m*t/(1+x);t2=m*(m-1)*t/(1+x)**2
   vals=[t,t1,t2]
   if k:
    g=lg-(1-x).log();vals=[t*g,t1*g+t/(1-x),t2*g+2*t1/(1-x)+t/(1-x)**2]
   for j in range(3):out[j]+=a*vals[j]
  return out
 return cs,norms

def analytic(F,K,p,U,normU):
 R=d.RA;k=F[1]
 S=(k.sqrt()*R).sinh()/k.sqrt() if k>0 else R
 G=((k.sqrt()*R).cosh()-1)/k if k>0 else R*R/2
 ball=None
 for i in range(1,1000):
  x=arb(i)/1000
  val=p*S+G*(F[0]+sum((F[j]*x**j for j in range(2,len(F))),arb(0)))
  contraction=G*sum((j*F[j]*x**(j-1) for j in range(2,len(F))),arb(0))
  if val<x and contraction<1:ball=x;break
 assert ball is not None,('no normalized complex source ball',str(p),[str(z) for z in F])
 kmax=peval(K,ball);fz=peval(F,ball)
 qz=peval([p*p]+[2*F[j-1]/j for j in range(1,7)],ball)
 u0,u1,u2=normU(ball)
 den=1-R*R*kmax/2;assert den>0
 h=(R*U[0]*p+R*R*(2*u0*fz+u1*(qz+2*p*p))/2)/den
 kd=peval([j*K[j] for j in range(1,len(K))],ball)*h+2*u0*kmax+2*u1*fz+u2*qz
 x=kmax.sqrt()*R;b=x.sinh()/x-1;assert b<1
 db=R*R*kd*(1+b)/(2*den)
 gap=arb('.002');rh=arb('.648');rm=d.RT
 bp=(b,db);bp1=d.ss(bp,1/gap);Db=d.ss(bp1,rh)
 DDb=d.ss(bp,rh/gap+2*rh*rh/gap**2)
 A=d.ss(d.sa(bp,Db),arb('.5'));DA=d.ss(d.sa(Db,DDb),arb('.5'))
 inv=(1/(1-b),db/(1-b)**2);i2=d.sm(inv,inv);i3=d.sm(i2,inv)
 C=(2/gap**2+1/(rm*gap),arb(0))
 C=d.sa(C,d.ss(d.sm(bp1,inv),1/gap))
 C=d.sa(C,d.ss((i2[0]-1,i2[1]),(rh/gap+2*rh**2/gap**2)/rm**2))
 C=d.sa(C,d.ss(d.sm(Db,i3),rh/gap/rm**2))
 C=d.sa(C,d.ss(d.sm(A,i2),2*rh/gap/rm**2))
 C=d.sa(C,d.ss(d.sm(DA,i2),1/rm**2))
 C=d.sa(C,d.ss(d.sm(d.sm(A,A),i2),1/rm**2))
 C=d.sa(C,d.ss(d.sm(d.sm(Db,A),i3),1/rm**2))
 C=d.sa(C,(kmax/2,kd/2))
 W0=((1-b)**(-arb('.5')),db/2*(1-b)**(-arb('1.5')))
 vals=[W0];W=W0
 for j in range(1,4):W=d.ss(d.sm(d.sm(W0,W0),d.sm(C,W)),arb(1)/j);vals.append(W)
 vals.append(d.sm(C,W))
 return dict(radius=ball,slack=ball-val,contraction=contraction,bmax=b,dbmax=db,kdelta=kd,Kmax=kmax,norms=vals,U_norms=[u0,u1,u2])

def transport(F,p,U,ab):
 N=d.N;K=[(j+1)*F[j+1] for j in range(5)]
 Z=[arb(0)]*(N+1);Z[1]=p
 for n in range(N-1):Z[n+2]=d.co(F,d.arb_poly(Z).truncate(n+1))[n]/((n+1)*(n+2))
 Z=d.arb_poly(Z);kz=d.co(K,Z);fz=d.co(F,Z)
 qz=d.co([p*p]+[2*F[j-1]/j for j in range(1,7)],Z)
 uz=d.co(U,Z);upz=d.co([(j+1)*U[j+1] for j in range(len(U)-1)],Z)
 uppz=d.co([(j+2)*(j+1)*U[j+2] for j in range(len(U)-2)],Z)
 src=2*d.mul(uz,fz)+d.mul(upz,qz+2*p*p)
 H=[arb(0)]*(N+1);H[1]=U[0]*p
 for n in range(N-1):H[n+2]=(sum((kz[j]*H[n-j] for j in range(n+1)),arb(0))+src[n])/((n+1)*(n+2))
 kd=d.mul(d.co([j*K[j] for j in range(1,len(K))],Z),d.arb_poly(H))+2*d.mul(uz,kz)+2*d.mul(upz,fz)+d.mul(uppz,qz)
 B=[arb(0)]*(N+1);B[1]=arb(1);C=[arb(0)]*(N+1)
 for n in range(N-1):
  B[n+2]=sum((kz[j]*B[n-j] for j in range(n+1)),arb(0))/((n+1)*(n+2))
  C[n+2]=sum((kz[j]*C[n-j]+kd[j]*B[n-j] for j in range(n+1)),arb(0))/((n+1)*(n+2))
 b=(d.arb_poly(B).right_shift(1)-1,d.arb_poly(C).right_shift(1));kp=(kz,kd)
 inv=d.po(b,arb(1));i2=d.pm(inv,inv);i3=d.pm(i2,inv);W0=d.po(b,arb('.5'));A=d.ps(d.pa(b,d.op(d.D,b)),arb('.5'))
 def op(w):
  flat=d.op(lambda a:d.arb_poly([(j+2)**2*a[j+2] for j in range(N-1)]),w)
  terms=[flat,d.pm(d.pm(d.op(d.df,b),inv),d.op(d.df,w)),
   d.pm(d.op(lambda a:d.sh(a,-2),(i2[0]-1,i2[1])),d.op(lambda a:d.D(d.D(a)),w)),
   d.pm(d.pm(d.op(lambda a:d.sh(d.D(a),-2),b),i3),d.op(d.D,w)),
   d.pm(d.pm(d.op(lambda a:d.sh(a,-2),d.ps(A,2)),i2),d.op(d.D,w)),
   d.pm(d.pm(d.op(lambda a:d.sh(d.D(a),-2),A),i2),w),
   d.pm(d.pm(d.op(lambda a:d.sh(a,-2),d.pm(A,A)),i2),w),
   d.pm(d.pm(d.pm(d.op(lambda a:d.sh(d.D(a),-2),b),A),i3),w),d.ps(d.pm(kp,w),arb('.5'))]
  out=(d.arb_poly([]),d.arb_poly([]))
  for x in terms:out=d.pa(out,x)
  return out
 W=W0;vals=[W]
 for j in range(1,4):
  W=d.pm(W0,d.op(lambda a:d.arb_poly([a[k]/(k+j) for k in range(N+1)]),d.pm(W0,op(W))))
  vals.append(W)
 vals.append(op(W));return vals

def table(vals,ab,g,defect=False):
 tab={key:[arb(0),arb(0)] for key in ['metric','mixed']}
 entries=[(3,vals[4],ab['norms'][4])] if defect else [(j,vals[j],ab['norms'][j]) for j in range(4)]
 for j,W,M in entries:
  for ix,which in enumerate(['inf','row']):
   def nm(s=0,p=0,k=0,shift=0):return d.norm(W[s],M[s],p,k,which,shift)
   def qw(k):return nm(1,0,k)+g*(nm(0,1,k)+nm(0,1,k,1)+2*nm(0,0,k)+nm(0,0,k,1))
   h=abs(j-1);fac=arb(1) if defect else d.T0**j
   tt={'metric':qw(0)+2*g*nm(k=1),'mixed':h*qw(0)+qw(1)+2*g*(h+1)*nm(k=1)+2*g*nm(k=2)}
   for key,val in tt.items():tab[key][ix]+=fac*val*(arb('1.14') if which=='row' else 1/(4*arb.pi()))
 return tab

def cutnorm(vals,ab,g):
 cg=1+g+ab['kdelta']/100+ab['dbmax']/100;out=[arb(0),arb(0)]
 for j in range(4):
  W=vals[j];M=ab['norms'][j]
  for ix,which in enumerate(['inf','row']):
   val=arb(0)
   for s,power in [(0,0),(0,1),(0,2),(1,0),(1,1)]:
    term=sum((W[s][n]*n**power*d.cutmoment(n,which) for n in range(d.KEEP+1)),arb(0))
    term+=M[s]/d.RT**(d.KEEP+1)*d.tailsum(power)*d.cutmoment(d.KEEP+1,which)
    val+=(cg if s==0 else 1)*term
   out[ix]+=d.T0**j*val
 for ix,which in enumerate(['inf','row']):out[ix]*=arb('1e8')*d.T0**-6*(arb('1.14') if which=='row' else 1/(4*arb.pi()))
 return out

def boundary(M,beta,gamma,eps,eta):
 V,t=sp.symbols('V t',nonnegative=True);q=sp.Rational
 f=q('1.06');G=q('.31')*V
 dx=[2*f*q('.30')*V**3,q('.9')*V**2];tx=[2*f*q('.31')*V**2,q('.61')*V]
 dtx=[2*(f*q('.9')*V**3+q('.3')*q('.31')*V**5),q('1.8')*V**2]
 sx=[2*f*3*V**4,3*V**3];dsx=[2*(f*3*V**5+q('.3')*3*V**7),3*V**4]
 tsx=[2*(f*3*V**4+q('.31')*3*V**6),3*V**3]
 dtsx=[2*(f*3*V**5+(q('.3')+q('.9'))*3*V**7)+2*q('.31')*3*V**7,3*V**4]
 ht={(1,0,0):q('.7'),(0,1,0):q('.7')*t*G,(2,0,0):q('1.2'),(1,1,0):q('1.2')*t*G,(0,2,0):q('.7')*t,
  (3,0,0):q('3.5'),(2,1,0):q('3.5')*t*G,(1,2,0):q('1.2')*t,(0,3,0):q('1.5')*t*t*G,
  (1,0,1):q('1.2')*G*G,(0,1,1):q('.7')*G,(2,0,1):q('3.5')*G*G,(1,1,1):q('1.2')*G,(0,2,1):q('.7')}
 def hd(I,j=None,k=None,nt=0):
  nf,ng=I
  for a in [j,k]:
   if a==0:nf+=1
   if a==1:ng+=1
  return ht[nf,ng,nt]
 def de(I):return sum(hd(I,j)*dx[j] for j in range(2))
 def dt(I):return sum(hd(I,j)*tx[j] for j in range(2))+t*hd(I,nt=1)
 def ddt(I):return sum(hd(I,j)*dtx[j]+t*hd(I,j,nt=1)*dx[j] for j in range(2))+sum(hd(I,j,k)*dx[k]*tx[j] for j in range(2) for k in range(2))
 val=V**4*(de((0,0))+ddt((0,0)))
 for j,I in enumerate([(1,0),(0,1)]):val+=hd(I)*dsx[j]+de(I)*sx[j]+dt(I)*dsx[j]+hd(I)*dtsx[j]+ddt(I)*sx[j]+de(I)*tsx[j]
 poly=sp.Poly(sp.expand(val),V,t);b1=arb(0);bi=arb(0)
 for (p,j),c in poly.terms():
  assert c>0 and p+2*j>=5
  cf=arb(int(c.p))/int(c.q)*eta**(arb(p-5)/2+j)
  b1+=cf*arb(2**(p+1)-1)/(p+1)*2*(arb(p+3)/2).gamma()
  bi+=cf*2**p*2*(arb(p)/2+1).gamma()
 b1*=beta.sqrt()/arb.pi()*M*gamma**3;bi*=arb('1.06')/arb.pi()*M*gamma**3
 c=arb(7)/8
 exited=220*beta.sqrt()/arb.pi()*(eta-c/(16*eta)).exp()*(arb(25)/(64*eta**4)+arb(95)/(16*c*eta**3)+8/(c*c*eta*eta))
 iexit=arb(2816)/arb.pi()*eta.exp()*eta**(-arb(7)/2)*(-3/(64*eta)).exp()
 b1+=eps*(exited*eta+arb('1e-100'));bi+=eps*(iexit*eta+arb('1e-100'))
 return (b1*bi).sqrt(),dict(polynomial=str(poly.as_expr()),L1_coefficient=enc(b1),Linf_coefficient=enc(bi))

def execute():
 ctx.prec=512;d.N=168;d.KEEP=144;d.RA=arb('.65');d.RT=arb('.64');d.RC=arb('.625');d.T0=arb('.0001')
 d.MOMENTS={};d.UPPERGAMMA={};cache={}
 def tails(p):
  if p not in cache:
   x=sp.symbols('x');f=x**(d.KEEP+1)/(1-x)
   for _ in range(p):f=x*sp.diff(f,x)
   a=sp.Rational((f/x**(d.KEEP+1)).subs(x,sp.Rational(125,128)))
   cache[p]=arb(int(a.p))/int(a.q)
  return cache[p]
 d.tailsum=tails
 old=result('current-localized-short-time-trial');p={k:dec(x) for k,x in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 L,gamma,beta,C,B2=[p[k] for k in ['L','gamma','beta','C','B2']];B0=L+gamma
 mp=result('current-global-mass-conversion')['parameters'];u,ep,tau=[dec(mp[k]) for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+ep,-2*ep,ep]
 # Convex Pseed on [-.4,1] and the exact g/D comparison imply qtilde<2.
 ratio=mx(arb(1),(1-L)/gamma)
 assert 6*beta+2*u-12*aa(ep)>0
 pm=beta-arb('.8')*beta+(3*beta+u)*arb('.16')+ep*arb('.4')**4
 assert ratio**2*mx(pm,arb(1))<2
 v=read(old['exact_trial']['Vpsi']);qb=read(old['exact_trial']['Q_over_gamma_squared'])
 qt=l.mul(l.mono(2),qb);wlo=gamma/(B0+arb('.4'));wzero=gamma/B0
 ev=[v]
 for _ in range(3):ev.append(E(ev[-1]))
 vv=[sup(a,wlo,arb(1),3) for a in ev]
 qe1=sup(E(qt),wlo,arb(1));qe2=sup(E(E(qt)),wlo,arb(1))
 jets=[vv[0],arb(2).sqrt()*vv[1],2*vv[2]+qe1*vv[1]/2,
  arb(2).sqrt()*(2*vv[3]+3*qe1*vv[2]/2+qe2*vv[1]/2)]
 alpha=3*mx(arb(1),(1-L)/gamma);localweight=(alpha*d.RC).exp();g=mx(jets[0],2*jets[1])*localweight
 worst={kind:{key:[arb(0),arb(0)] for key in ['metric','mixed']} for kind in ['P','F']}
 wcap=(gamma/(B0-arb('.995')))**3
 eps=mx(*(jets[j]*wcap/arb([1,10,100,10000][j]) for j in range(4)))
 for length in [arb('.5'),arb(2)]:
  a,b,c=[jets[j]*wcap for j in range(3)];s=b+length*c
  ext=[a+length*b+length*length*c/2,s,c+arb('1.875')/length*s,
       2*arb('1.875')/length*c+6/length**2*s]
  eps=mx(eps,*(ext[j]/arb([1,10,100,10000][j]) for j in range(4)))
 for i in range(64):
  lo=arb(i)/64;hi=arb(i+1)/64 if i<63 else L
  z=(lo+hi)/2+arb(0,(hi-lo)/2);dm=B0-lo
  P=beta+2*beta*z+(3*beta+u)*z*z+ep*z**4
  ps=ratio*up(P).sqrt()
  F=[ratio*aa(z*((6*beta+2*u)*z-u+ep*(3*z**3-2*z*z)))+3*tau/gamma*aa(z*z)/2]
  for j in range(1,6):F.append(dm**(j-1)*aa(sum((arb(k)*Q[k]*comb(k-1,j)*z**(k-1-j)/2 for k in range(j+1,7)),arb(0))))
  K=[(j+1)*F[j+1] for j in range(5)]
  a=gamma/(B0-lo);b=gamma/(B0-hi) if i<63 else arb(1)
  U,normU=Udata(v,arb(a.lower()),arb(b.upper()))
  ab=analytic(F,K,ps,U,normU);vals=transport(F,ps,U,ab)
  tabs={'P':table(vals,ab,g),'F':table(vals,ab,g,True)};cut=cutnorm(vals,ab,g)
  for kind in ['P','F']:
   for key in ['metric','mixed']:
    for j in range(2):
     tabs[kind][key][j]=localweight*(tabs[kind][key][j]+cut[j]*(d.T0**3 if kind=='P' else 1))
     worst[kind][key][j]=mx(worst[kind][key][j],tabs[kind][key][j])
  CELLS.append(dict(index=i,z_interval=[str(lo),str(hi)],source_ball=enc(ab['radius']),source_ball_slack=enc(ab['slack']),
   source_contraction=enc(ab['contraction']),analytic_metric_curvature=enc(ab['kdelta']),analytic_metric_Jacobi=enc(ab['dbmax']),
   source_tables={kind:{key:[enc(x) for x in pair] for key,pair in tab.items()} for kind,tab in tabs.items()}))
  if i%4==3:print(json.dumps({'completed_source_cells':i+1,'total':64}),flush=True)
 raw=json.loads((ROOT/'evidence/current-disk-heat-propagation.json').read_bytes())['runs'][-1]
 rows={n['kind']:n for n in raw['numerics']};cap=rows['conditional_complete_cap_kernel_table']['table']
 tables=rows['full_uncut_transport_norms']['tables'];cut=rows['complete_actual_domain_propagation_budget']['cutoff_jet_norm']
 P0={};F0={}
 for key in ['base','time']:
  a=pair(tables['P'][key]);b=pair(cap[key]);f=pair(tables['F'][key])
  for j,w in enumerate(['inf','row']):a[j]+=dec(cut[w])*arb('.001')**3;f[j]+=dec(cut[w])+arb('1e-25')+arb('1e-100')
  b[0]+=arb('1e-100');P0[key]=[mx(2*b[0],a[0]),b[1]+a[1]];F0[key]=f
 for key in ['metric','mixed']:
  c=pair(cap[key]);scaled=[2*c[0]*eps/wzero**3,c[1]*eps/wzero**3]
  worst['P'][key]=[mx(worst['P'][key][0],scaled[0]),worst['P'][key][1]+scaled[1]]
  worst['F'][key]=[x+eps/wzero**3*(arb('1e-25')+arb('1e-100')) for x in worst['F'][key]]
 keys=['A_delta','B_delta','A_deltaD','B_deltaD','a_delta','b_delta','a_deltaD','b_deltaD']
 kvals=worst['P']['metric']+worst['P']['mixed']+worst['F']['metric']+worst['F']['mixed']
 wc=(4*alpha).exp();PW={k:[wc*x for x in p] for k,p in P0.items()};FW={k:[wc*x for x in p] for k,p in F0.items()}
 Pd,PdT,Fd,FdT=worst['P']['metric'],worst['P']['mixed'],worst['F']['metric'],worst['F']['mixed']
 e0=cv(PW['base'],FW['base']);ed=cv(Pd,FW['base'])+cv(PW['base'],Fd)
 et=e0+cv(PW['time'],FW['base'])+cv(PW['base'],FW['time'])
 edt=ed+cv(PdT,FW['base'])+cv(PW['time'],Fd)+cv(Pd,FW['time'])+cv(PW['base'],FdT)
 x0=d.T0**4/140*FW['base'][1];xt=d.T0**4/140*(FW['base'][1]+FW['time'][1]);xd=d.T0**4/140*Fd[1];xdt=d.T0**4/140*(Fd[1]+FdT[1])
 assert x0<1
 ret=4*(edt/(1-x0)+ed*xt/(1-x0)**2+et*xd/(1-x0)**2+e0*(xdt/(1-x0)**2+2*xd*xt/(1-x0)**3))
 capbase,capT=pair(cap['base']),pair(cap['time'])
 capd=[eps*x for x in pair(cap['metric'])];capdT=[eps*x for x in pair(cap['mixed'])]
 fscale=wzero**3*(alpha*(2+d.RC)).exp();sfd=[fscale*x for x in Fd];sfdT=[fscale*x for x in FdT]
 smooth=4*(cv(capd,F0['base'])+cv(capbase,sfd)+cv(capdT,F0['base'])+cv(capT,sfd)+cv(capd,F0['time'])+cv(capbase,sfdT))
 zabs=arb('.4');dmin=B0-zabs
 qn=[sum((aa(arb(comb(j,k))*Q[j])*zabs**(j-k) for j in range(k,7)),arb(0)) for k in range(3)]
 qn[2]*=2
 mv=arb('1.5')*((qn[1]*zabs**5+5*qn[0]*zabs**4)/dmin**4+4*qn[0]*zabs**5/dmin**5)+zabs**5*(zabs*qn[2]+4*qn[1])/dmin**4
 md=(arb('1.72')+arb('1.62')*arb('.75'))*mv
 M=arb(0);xx=sp.symbols('x')
 for n0,major in [(4,mv),(4,md),(5,md)]:
  for j in range(3):
   expr=sp.diff(xx**n0/(1-xx),xx,j)/xx**(n0-j)
   a=sp.Rational(expr.subs(xx,sp.Rational(2,3)))
   M=mx(M,major/arb('.75')**n0*arb(int(a.p))/int(a.q))
 b7,bdetail=boundary(M,beta,gamma,eps,d.T0)
 plus=result('current-positive-time-quantum-operator')['constants'];CL=dec(plus['CL_t']);kp=dec(plus['total_operator_difference_upper'])
 Elower=dec(old['constants']['trial_Et_lower']);vnorm=dec(old['constants']['original_v_L2'])
 pref=4*C*CL/(B2*beta)/Elower;Wnorm=(gamma*(1-wzero**5)/5).sqrt()
 tauprop=pref*Wnorm*ret*d.T0**2/2;tausmooth=pref*smooth*d.T0**2/2
 taub=pref*b7*arb(4)/3*d.T0**(arb(3)/4)
 tau6=dec(old['constants']['attained_leading_column_upper'])*d.T0/arb('.001')
 layer=result('current-localized-boundary-layer-test')['constants'];tau5=dec(layer['explicit_layer_relative_column_upper'])*(d.T0/arb('.001'))**(arb(3)/4)
 D=lambda f:l.mul(l.mono(2),l.diff(f));Delta=lambda f:D(l.mul(qb,D(f)))
 R=l.scale(D(D(qb)),-1);dv=Delta(v)
 p4=l.add(l.scale(l.mul(l.mul(R,R),v),4),l.scale(l.mul(R,dv),8),l.scale(Delta(l.mul(R,v)),8),l.scale(Delta(dv),8))
 p4norm=l.l2(p4,wzero,gamma)
 eold=rows['complete_actual_domain_propagation_budget'];mh=4*P0['base'][0]+(dec(eold['first_actual']['base'])+dec(eold['higher_Volterra']))*arb('.001')**4/4
 window=4*C*CL/(B2*beta)*vnorm/Elower*((64*mh+2/arb.pi())*(d.T0**-2-arb('.001')**-2)/2+p4norm/vnorm/(480*arb.pi())*(arb('.001')/d.T0).log())
 total=tauprop+tausmooth+taub+tau6+2*tau5+window
 target=dec(old['constants']['conditional_relative_column_screen'])
 Ghat=dec(result('current-first-quantum-solve-attempt')['constants']['proved_source_dual_upper'])
 Kfirst=dec(result('current-first-inverse-source-screen')['constants']['first_residual_to_original_bolt_multiplier'])
 const=dict(eta=d.T0,global_weighted_v=jets[0],global_weighted_gradient=jets[1],fixed_endpoint_multiplier=g,
  cap_epsilon=eps,cap_epsilon_over_W0=eps/wzero**3,local_boundary_Cauchy_coefficient=M,
  boundary_L2_t7over4_coefficient=b7,complete_return_coefficient=ret,Volterra_base_ratio=x0,
  tau_propagation=tauprop,tau_smooth_cap=tausmooth,tau_boundary=taub,tau_leading=tau6,
  tau_explicit_layer_overcount=2*tau5,tau_full_extra_window=window,full_relative_column_upper=total,
  conditional_column_screen=target,screen_fraction=total/target,
  original_first_error_envelope=Kfirst*Ghat*total/(1-kp),P4_original_L2=p4norm,window_heat_diagonal_coefficient=mh)
 assert all(x.is_finite() for x in const.values())
 return dict(precision_bits=ctx.prec,retained_coefficients=d.KEEP,computed_coefficients=d.N,
  constants={k:enc(v) for k,v in const.items()},actual_kernel_inputs={k:enc(x) for k,x in zip(keys,kvals)},
  local_boundary_argument=bdetail,column_screen_pass=bool(total<target),
  actual_Ritz_coefficient=None,complete_first_residual=None,signed_Gt=None,
  scope=dict(actual_paired_trial=True,full_angular_and_thermal_sums=True,full_analytic_tails=True,
   common_regulator_limit_for_this_column=True,canonical_massive_column_only=True,
   actual_first_correction=False,full_coupled_admission=False,nonlinear_existence_matching_validity_stability=False,goal_complete=False))

if __name__=='__main__':
 out=ROOT/'evidence'/(STEM+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result']=execute();after={p:digest(ROOT/p) for p in INPUTS};assert before==after
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(status='complete_column_arithmetic_pending_argument_review',input_sha256_after=after,protected_inputs_unchanged=len(protected))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec.update(finished_utc=datetime.now(timezone.utc).isoformat(),source_cells=CELLS)
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),completed_cells=len(CELLS),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if rec['status']=='failed':raise SystemExit(1)
