"""Actual-source, all-coefficient lower test of the full first residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
from flint import acb
import current_localized_short_time_trial as lp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-first-residual-range-test'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_localized_first_residual_range_test.py',
 'checks/current_quartic_ritz_rank_test.py','checks/current_localized_short_time_trial.py',
 'notes/current-finite-local-action-assembly.md',
 'notes/current-canonical-heat-restoration-assembly.md',
 'notes/current-quartic-ritz-rank-test.md',
 'notes/current-scalar-green-source-residual.md',
 'notes/current-positive-time-quantum-operator.md',
 'notes/current-localized-quantum-source-test.md',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-first-inverse-source-screen.json',
 'evidence/current-first-quantum-solve-attempt.json',
 'evidence/current-scalar-green-source-residual.json',
 'evidence/current-positive-time-quantum-operator.json',
 'evidence/current-localized-short-time-trial.json',
 'evidence/current-localized-covariant-jet-test.json',
 'evidence/current-localized-quantum-source-test.json',
 'evidence/current-localized-first-residual-range-test-attempt1.json',
 'evidence/current-localized-first-residual-range-test-attempt2.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(n):
 d=json.loads((ROOT/'evidence'/(n+'.json')).read_bytes())
 r=d.get('runs',[d])[-1];assert r['status']!='failed'
 return r['result']
def constants(n):return {k:dec(v) for k,v in result(n)['constants'].items()}
def upper(x):return arb(arb(x).upper())
def mag(x):return upper(abs(x))
def rad(x):return arb(0,mag(x))
def read(items):return {(x['power'],x['log_power']):dec(x['coefficient']) for x in items}

def pa(*ps):
 r=[arb(0)]*max(map(len,ps))
 for p in ps:
  for i,v in enumerate(p):r[i]+=v
 return r
def ps(p,c):return [v*c for v in p]
def pm(p,q):
 r=[arb(0)]*(len(p)+len(q)-1)
 for i,v in enumerate(p):
  for j,w in enumerate(q):r[i+j]+=v*w
 return r
def pp(p,n):
 r=[arb(1)]
 for _ in range(n):r=pm(r,p)
 return r
def pd(p):return [i*p[i] for i in range(1,len(p))] or [arb(0)]
def pv(p,x):
 r=arb(0)
 for v in reversed(p):r=r*x+v
 return r
def pint(p):return [arb(0)]+[v/(i+1) for i,v in enumerate(p)]
def jets(p,x,n=3):
 r=[]
 for _ in range(n+1):r.append(pv(p,x));p=pd(p)
 return r

class Jet:
 def __init__(self,v=0,i=0,j=0,ij=0):self.v,self.i,self.j,self.ij=map(arb,[v,i,j,ij])
 @staticmethod
 def cast(x):return x if isinstance(x,Jet) else Jet(x)
 def __add__(self,o):
  o=self.cast(o);return Jet(self.v+o.v,self.i+o.i,self.j+o.j,self.ij+o.ij)
 __radd__=__add__
 def __neg__(self):return Jet(-self.v,-self.i,-self.j,-self.ij)
 def __sub__(self,o):return self+-self.cast(o)
 def __rsub__(self,o):return self.cast(o)+-self
 def __mul__(self,o):
  o=self.cast(o)
  return Jet(self.v*o.v,self.i*o.v+self.v*o.i,self.j*o.v+self.v*o.j,
   self.ij*o.v+self.i*o.j+self.j*o.i+self.v*o.ij)
 __rmul__=__mul__
 def __truediv__(self,o):return self*(1/arb(o))
 def __pow__(self,n):
  r=Jet(1)
  for _ in range(n):r=r*self
  return r
def ex(i,j,k):return Jet(1,k*i,k*j,k*k*i*j)

def local_density(p,z,qj,di,dj):
 gamma,C,B2=p['gamma'],p['C'],p['B2']
 q,qx,qxx=qj[:3];ni,ui=di;nj,uj=dj
 delta=lambda f:q*f[2]+qx*f[1]
 phi=Jet(-z.log(),ni[0],nj[0]);L=Jet(-qxx,-2*delta(ui),-2*delta(uj))
 M=Jet(gamma*qx/z+gamma**2*q/z**2,delta(ni),delta(nj))
 V=Jet(gamma**2*q/z**2,2*gamma*q*ni[1]/z,2*gamma*q*nj[1]/z,2*q*ni[1]*nj[1])
 D=L-6*M-6*V
 n2,n4=ex(ni[0],nj[0],2),ex(ni[0],nj[0],4)
 u2,um2=ex(ui[0],uj[0],2),ex(ui[0],uj[0],-2)
 le=n2*D/(2*z**2)+n2*u2/(C*z**2)+3*arb(2000)**-10*B2/C*n4*u2/z**4
 lem=-arb.pi()*arb(2000)**16/(C*B2)*u2
 lr=8*arb.pi()*arb('1e22')/(C*B2)*(4*u2+4*C*D+C*C*um2*D**2)
 wb=phi*L**2/60+L*M/30-M**2/10+L*V/30-7*M*V/90-7*V**2/180
 lw=(C*C*um2*wb-C*phi*L/18+(arb(1)/15-arb(2000)**8/3)*phi*u2-C*V/18)/(2*arb.pi()*C*B2)
 return (le+lem+lr+lw).ij

def qreal(k,w):
 if k>100:
  # Complete positive-integral bound, avoiding the large-degree library NaN.
  # Q_k(w)<=Q_0(w)*(w+sqrt(w^2-1))^-k for k>=0, w>1.
  return rad(((w+1)/(w-1)).log()/2*(-k*(w+(w*w-1).sqrt()).log()).exp())
 v=acb(w).legendre_q(k,type=3)
 assert v.is_finite() and v.imag.contains(0),(k,w,v)
 return v.real

def qmoment(f,k,a,b):
 lam=a*k*(k+1);out=[arb(0)]*len(f)
 for n in range(len(f)-1,-1,-1):
  den=a*n*(n+1)-lam;assert not den.contains(0)
  out[n]=(f[n]-(2*(n+1)**2*out[n+1] if n+1<len(f) else 0))/den
 for n in range(len(f)):
  rem=(a*n*(n+1)-lam)*out[n]+(2*(n+1)**2*out[n+1] if n+1<len(f) else 0)-f[n]
  assert rem.contains(0)
 w=1+a*b;q=qreal(k,w);flux=k*(w*q-qreal(k-1,w))
 return q*(2*b+a*b*b)*pv(pd(out),b)-pv(out,b)*flux-out[0]

def execute():
 ctx.prec=512
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 L,gamma,a,beta,C,B2=[p[k] for k in ['L','gamma','a','beta','C','B2']]
 scalar=result('current-scalar-green-source-residual');sc=constants('current-scalar-green-source-residual')
 fs=constants('current-first-inverse-source-screen');pl=constants('current-positive-time-quantum-operator')
 col=constants('current-localized-quantum-source-test')
 trial=result('current-localized-short-time-trial');tc=constants('current-localized-short-time-trial')
 At=dec(result('current-global-green-trial-enclosure')['constants']['A_t'])
 y=[(dec(v['degree']),dec(v['coefficient'])/At) for v in scalar['profiles']['y0']]
 assert At>0 and beta>arb('.16') and L<=1
 assert C>1 and arb(2000)**-10*B2/C<1 and arb.pi()*arb(2000)**16/(C*B2)<1
 assert 8*arb.pi()*arb('1e22')*C/B2<1
 assert C/(480*arb.pi()*B2)<arb('.01')
 assert qreal(arb(0),arb(2)).overlaps(arb(3).log()/2)
 mp=result('current-global-mass-conversion')['parameters']
 u,ep,tau=[dec(mp[k]) for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+ep,-2*ep,ep]
 q=[arb(0),arb(2),a]
 for j in range(3,7):q.append((-1)**j*sum((comb(i,j)*Q[i]*L**(i-j) for i in range(j,7)),arb(0))*gamma**(j-2))
 z=[L,-gamma];z4=pp(z,4)
 # Original parameter guards, not a midpoint finite geometry.
 assert sum((abs(c)*L**max(0,i-1)*i for i,c in enumerate(Q)),arb(0))<10
 assert sum((abs(c)*L**max(0,i-2)*i*(i-1) for i,c in enumerate(Q)),arb(0))<30
 esy=fs['scalar_source_uniform_upper']*fs['actual_trial_Hmu_norm_upper']+sc['full_scalar_residual_Es_dual_upper']
 S2=4*C*C/arb(2000)**4
 # S2 denominator q_ang=2000^4 (not q_ang squared).
 F_eff=rad(100);c0=2*(604*arb(2000)**8+(1+2*arb('1e-8'))*(arb(2000)**8/6-arb(1)/30))/(arb.pi()*C*B2)
 assert c0<1
 # Broad a2 ball is safe here: it is a known interval, not a zeroed contact.
 b2=C*F_eff/(480*arb.pi()*B2)
 epsbar=C/(60*arb.pi()*B2);Cv=1+11*epsbar*beta/3
 # Exact psi global weighted energy comparison, using its z^4 factors.
 qb=read(trial['exact_trial']['Q_over_gamma_squared']);zw=lp.add(lp.mono(0,L+gamma),lp.mono(-1,-gamma))
 w0=dec(trial['exact_trial']['w0']);dw=lambda f:lp.mul(lp.mono(2),lp.diff(f))
 dz=lambda f:lp.scale(dw(f),1/gamma)
 Qw=lp.scale(qb,gamma**2);Qpw=lp.scale(dw(qb),gamma)
 If=lp.primitive_from(lp.mul(lp.add(lp.mul(zw,dw(dw(qb))),lp.scale(Qpw,-4)),lp.mono(2)),w0)
 Ab=lp.add(lp.mul(lp.mul(zw,dw(qb)),lp.mono(4)),lp.scale(lp.mul(qb,lp.mono(4)),5*gamma),lp.scale(lp.mul(lp.mul(qb,zw),lp.mono(5)),4))
 PsiBar=lp.add(Ab,lp.scale(If,-1))
 NuBar=lp.add(lp.scale(PsiBar,arb(1)/2),lp.scale(lp.mul(qb,lp.mono(4)),gamma))
 VBar=lp.add(lp.scale(Ab,-arb(3)/2),If)
 def intw(f):return gamma*sum((v*lp.moment(n-2,k,w0) for (n,k),v in f.items()),arb(0))
 Jpsi=arb(0)
 for f in [NuBar,VBar]:
  g0=lp.mul(lp.power(zw,2),f)
  g1=lp.add(lp.scale(lp.mul(zw,f),4),lp.mul(lp.power(zw,2),dz(f)))
  g2=lp.add(lp.mul(Qpw,lp.add(lp.scale(lp.mul(lp.power(zw,2),f),4),lp.mul(lp.power(zw,3),dz(f)))),
   lp.mul(Qw,lp.add(lp.scale(lp.mul(zw,f),12),lp.scale(lp.mul(lp.power(zw,2),dz(f)),8),lp.mul(lp.power(zw,3),dz(dz(f))))))
  for integrand in [lp.mul(g0,g0),lp.mul(Qw,lp.mul(g1,g1)),lp.mul(g2,g2)]:
   v=intw(integrand);assert v>=0,v;Jpsi+=v
 h3=constants('current-localized-covariant-jet-test')['covariant_3_L2']
 egeo_psi=4*S2/(525*arb.pi()*B2*beta)*h3*h3
 Epsi=upper(((1+arb('1e-70'))/(1-arb('.0022'))*(arb('1e8')*Jpsi+egeo_psi)).sqrt())
 assert Epsi<arb('1e-5'),Epsi
 # Canonical fields and x derivatives, exact finite Laurent/log expressions.
 psif=[read(trial['exact_trial'][key]) for key in ['Upsi','Vpsi']]
 dx=lambda f:lp.scale(dw(f),-1)
 psij=[]
 for f in psif:
  row=[]
  for j in range(4):row.append(f);f=dx(f)
  psij.append(row)
 coeff_tau=[u+3*beta+ep,-5*ep,arb(9)*ep/2]
 tau_base=pa([coeff_tau[0]],ps(z,coeff_tau[1]),ps(pp(z,2),coeff_tau[2]))
 out=[]
 for bint in [1,4]:
  b=arb(bint);assert b*gamma<sc['cutoff_radius']
  chi=[arb(comb(8,j))*(-1)**j/b**j for j in range(9)]
  qc=pm(q,chi)
  Ip=pint(pm(pa(pm(z,pd(pd(q))),ps(pd(q),4*gamma)),chi))
  I=pa([pv(Ip,b)],ps(Ip,-1))
  Abar=pa(ps(qc,5*gamma),ps(pm(z,pd(qc)),-1))
  Xibar=pa(Abar,ps(I,-1));xi=pm(z4,Xibar)
  nu=pa(ps(xi,arb(1)/2),ps(pm(z4,qc),gamma))
  vv=pm(z4,pa(ps(Abar,-arb(3)/2),I))
  # Endpoint and York polynomial identities, before norm arithmetic.
  for f in [xi,nu,vv]:
   for j in range(6):assert pv(f,b).contains(0);f=pd(f)
  f0=pm(tau_base,Xibar)
  core=pa(ps(f0,40*gamma**2),ps(pm(z,pd(f0)),-14*gamma),pm(pp(z,2),pd(pd(f0))))
  drift=pa(ps(pm(z,f0),-8*gamma),pm(pp(z,2),pd(f0)))
  f=ps(pm(z4,pa(ps(f0,Cv),ps(pa(pm(q,core),pm(pd(q),drift)),-p['d']))),4/beta)
  source=arb(0);source_parts=[]
  Q0int=((b+2/a)*(b+2/a).log()-b*b.log()-(2/a)*(2/a).log())/2
  for index,(k,c) in enumerate(y):
   if index==0:
    amount=mag(c)*sum((mag(v)*b**j for j,v in enumerate(f)),arb(0))*Q0int
    v=rad(amount)
   else:v=c*qmoment(f,k,a,b)
   source+=v;source_parts.append(enc(v))
  # Complete cell bounds. Precompute polynomial derivative arrays.
  fj=[[f,pd(f),pd(pd(f)),pd(pd(pd(f)))] for f in [nu,vv]]
  qj=[q,pd(q),pd(pd(q)),pd(pd(pd(q)))];xij=[xi,pd(xi)]
  ncell=2048*bint;step=b/ncell
  cross=arb(0);vac_abs=arb(0);Jtest=arb(0);geoint=arb(0);xiH2=arb(0)
  for n in range(ncell):
   x=arb((n+arb(1)/2)*step,step/2);zz=L-gamma*x
   qv=[pv(f,x) for f in qj]
   dirs=[[pv(f,x) for f in row] for row in fj]
   w=1/(1+x);pdirs=[[lp.val(f,w) for f in row] for row in psij]
   density=-local_density(p,zz,qv,pdirs,dirs)/beta
   vi,vj=pdirs[1],dirs[1];R=-qv[2]
   dvi=qv[0]*vi[2]+qv[1]*vi[1];dvj=qv[0]*vj[2]+qv[1]*vj[1]
   H4=4*R*R*vi[0]*vj[0]+8*R*(vi[0]*dvj+vj[0]*dvi)+8*dvi*dvj
   density-=b2*H4/beta+arb(2000)**4/(6*arb.pi()*B2*beta)*qv[0]*vi[1]*vj[1]
   density+=pl['a0_coefficient']*vi[0]*vj[0]+pl['a2_coefficient']*H4
   cross+=density*step;vac_abs+=mag(vi[0]*vj[0])*step*c0/beta
   for fsj in dirs:
    dl=qv[0]*fsj[2]+qv[1]*fsj[1]
    Jtest+=(mag(fsj[0]/zz**2)**2+mag(qv[0])*mag(fsj[1]/zz**2)**2+mag(dl/zz)**2)*step*gamma
   A=qv[0]*dirs[1][3];D=(3*qv[1]*dirs[1][2]+qv[2]*dirs[1][1])/2
   geoint+=mag(qv[0])*(arb(5)/8*mag(A)**2+arb(3)/2*mag(A*D)+mag(D)**2)*step*gamma
   xiH2+=mag(pv(xi,x)/zz**2)**2*step*gamma
  Etest=upper(((1+arb('1e-70'))/(1-arb('.0022'))*(arb('1e8')*Jtest+2*S2/(525*arb.pi()*B2*beta)*geoint)).sqrt())
  mixed_error=upper(sc['complete_mixed_remainder_constant']*esy*xiH2.sqrt()/(beta*At*gamma))
  short_error=upper(col['full_relative_column_upper']*tc['trial_Et_lower']*Etest/gamma)
  plus_error=upper(pl['positive_time_operator_norm_upper']*Epsi*Etest/gamma)
  transport_error=upper(arb('1e-70')*Epsi*Etest/gamma)
  Dfull=cross+rad(vac_abs+short_error+plus_error+transport_error)
  Gfull=source+rad(mixed_error)
  item=dict(b=bint,cells=ncell,source_signed=enc(Gfull),source_main=enc(source),source_parts=source_parts,
   mixed_remainder_error=enc(mixed_error),response_local=enc(cross),response_signed=enc(Dfull),
   vacuum_error=enc(vac_abs),short_error=enc(short_error),positive_time_error=enc(plus_error),
   retained_transport_error=enc(transport_error),
   test_energy_upper=enc(Etest),weighted_J=enc(Jtest),full_geodesic_energy_upper=enc(2*S2/(525*arb.pi()*B2*beta)*geoint))
  out.append((Gfull,Dfull,Etest,item))
  print(json.dumps({'test_b':bint,'G':str(Gfull),'D':str(Dfull),'energy_upper':str(Etest)}),flush=True)
 G1,D1,E1,_=out[0];G2,D2,E2,_=out[1]
 determinant=D2*G1-D1*G2
 denom=mag(D2)*E1+mag(D1)*E2
 lower=(arb(abs(determinant).lower())*arb(At.lower())*arb(gamma.lower())/upper(denom)) if not determinant.contains(0) else arb(0)
 lower=arb(lower.lower());screen=arb('1e33');bolt=arb('1.667e-9')
 floor=fs['first_residual_to_original_bolt_multiplier']*lower
 status='actual_source_one_column_residual_screen_rejected' if lower>screen else 'actual_source_range_test_inconclusive'
 # An explicit actual candidate exists even though (7) rejects every coefficient.
 center=(G1/D1).mid() if not D1.contains(0) else arb(0)
 cs=dict(source_amplitude=At,physical_measure_gamma=gamma,psi_weighted_J=Jpsi,
  psi_full_energy_upper=Epsi,source_scalar_energy_upper=esy,
  complete_signed_determinant=determinant,complete_first_residual_lower=lower,
  residual_lower_over_screen=lower/screen,first_output_estimator_floor=floor,
  estimator_floor_over_entire_bolt_allowance=floor/bolt,
  explicit_trial_coefficient_over_At=arb(center))
 return status,dict(precision_bits=ctx.prec,constants={k:enc(v) for k,v in cs.items()},tests=[x[3] for x in out],
  actual_source_not_arbitrary_input=True,complete_quantum_errors_retained=True,
  complete_residual_upper=None,signed_Gt=None,original_physical_error_lower=None,
  original_admission_closed=False,construction_obstruction=False,goal_complete=False,
  alternative_frozen_contract={k:enc(v) for k,v in constants('current-first-quantum-solve-attempt').items()
   if k in ['exact_application_first_residual_error','conditional_common_stage_output_error_ceiling',
            'first_inverse_contraction_upper','exact_application_analytic_physical_debit']},
  alternative_applications_executed=False)

if __name__=='__main__':
 path=ROOT/'evidence'/(STEM+'.json');assert not path.exists(),'Preserve immutable receipt.'
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
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 path.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
