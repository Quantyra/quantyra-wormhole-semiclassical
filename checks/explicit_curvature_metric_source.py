"""Exact charged metric source; independent 4D curvature and moment fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from itertools import product
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='explicit-curvature-metric-source'
OWN=['checks/explicit_curvature_metric_source.py','notes/'+STEM+'.md',
     'evidence/curvature-metric-source-review.json']
PREVIOUS='physical-curvature-source-and-boundary-range'
ROWS=[]; NUMERICS=[]
PROTOCOL={'scope':'Exact first-order metric source and its adjoint moment; physical constrained value unresolved',
 'mpmath_precision':55,'quadrature_relative_tolerance':'1e-38',
 'curvature_arithmetic':'exact rational first-variation jets; full four-dimensional contractions',
 'curvature_fixture_count':4,'moment_fixture_count':3,
 'moment_interval':['1.3','2.2'],'moment_seed':{'B':'1','lambda':'2','d':'0.01','horizon':'1.1'},
 'test_weights_are_horizon_solutions':False,'fixtures_are_physical_solutions':False,
 'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}; ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 val=S.simplify(S.sympify(value).doit());check(name,val==0,val)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==123 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('source_scope',len(review['primary_sources'])==1 and not review['primary_sources'][0]['Gauss_Bonnet_action_adopted']
       and not review['unconstrained_profiles_are_physical_solutions'])

def radial(F,r):
 R=-S.diff(F,r,2)-4*S.diff(F,r)/r-2*(F-1)/r**2
 W=-S.diff(F,r,2)+2*S.diff(F,r)/r+2*(1-F)/r**2
 j=-S.diff(F,r,2)/4+(F-1)/(2*r*r)
 c=7*W/18+11*j/3
 d=7*W/(3*r)+44*j/(3*r)+S.diff(R,r)
 e=44*j/(3*r*r)+S.diff(R,r,2)+2*S.diff(R,r)/r
 return R,W,j,c,d,e

def predictions(F,a,b,t,r):
 R,W,j,c,d,e=radial(F,r)
 dH=S.diff(b,r)+F*S.diff(a,r)
 dR2=-S.diff(b,r,2)-3*S.diff(F,r)*S.diff(a,r)-2*F*S.diff(a,r,2)-S.diff(b,t,2)/F**2
 psi=dR2-2*b/r**2-4*dH/r
 dC=2*W*(dR2+2*dH/r-2*b/r**2)/3
 dJ=2*j*(dR2+2*b/r**2)
 KR=b*S.diff(R,r,2)+(dH+2*b/r)*S.diff(R,r)
 return R,W,j,psi,dC,dJ,KR,c*psi+d*dH+e*b

def exact_checks():
 r,t=S.symbols('r t',positive=True);F=S.Function('F')(r)
 a=S.Function('a')(t,r);b=S.Function('b')(t,r)
 R,W,j,psi,dC,dJ,KR,mix=predictions(F,a,b,t,r)
 zero('exact_metric_coupling_rearrangement',KR+7*dC/12+11*dJ/6-mix)
 zero('scalar_contains_required_time_jet',S.diff(psi,S.diff(b,t,2))+1/F**2)
 y=S.symbols('y',positive=True);ay,by,beta,c,d,e,E,Z=S.symbols('ay by beta c d e E Z')
 br=-by+S.diff(F,r)*beta
 A1=r*F*d;A0=r*F*(d*S.diff(F,r)+e*F)
 zero('relative_metric_source_transform',r*F*(d*(br-ay)+e*F*beta)-(-A1*(ay+by)+A0*beta))
 uu,ps,gsrc,fA=S.symbols('uu ps gsrc fA')
 zero('diagonal_shift_and_source_sign',(-r*F*c*ps+gsrc-r*F*fA)-(gsrc-r*F*fA-F*c*uu).subs(uu,r*ps))
 lam,m,Q,c3,c4,B,ds,z,x,mu,q=S.symbols('lam m Q c3 c4 B ds z x mu q',positive=True)
 Fr=lam*r*r+1-2*m/r+Q/r**2+c3/r**3+c4/r**4
 R,W,j,c,d,e=radial(Fr,r)
 zero('seed_Weyl_scalar',W-(12*m/r**3-12*Q/r**4-20*c3/r**5-30*c4/r**6))
 zero('seed_tracefree_Ricci_scalar',j-(-Q/r**4-5*c3/(2*r**5)-9*c4/(2*r**6)))
 zero('seed_scalar_derivative',S.diff(R,r)-(10*c3/r**6+36*c4/r**7))
 sub={r:B/z,lam:x/B**2,m:mu*B,Q:q*B**2,c3:-2*ds*B**3,c4:ds*B**4}
 f=x+z*z-2*mu*z**3+q*z**4-2*ds*z**5+ds*z**6
 D=28*mu-S.Rational(128,3)*q*z+S.Rational(440,3)*ds*z*z-100*ds*z**3
 H=-S.Rational(44,3)*q+S.Rational(460,3)*ds*z-246*ds*z*z
 cp=(S.Rational(14,3)*mu*z**3-S.Rational(25,3)*q*z**4+S.Rational(305,9)*ds*z**5-S.Rational(169,6)*ds*z**6)/B**2
 zero('seed_f_polynomial',z*z*Fr.subs(sub,simultaneous=True)-f)
 zero('seed_c_polynomial',c.subs(sub,simultaneous=True)-cp)
 zero('seed_d_polynomial',d.subs(sub,simultaneous=True)-z**4*D/B**3)
 zero('seed_e_polynomial',e.subs(sub,simultaneous=True)-z**6*H/B**4)
 A1=r*Fr*d;A0=r*Fr*(d*S.diff(Fr,r)+e*Fr)
 P1=z*f*D/B**2;P0=f*(D*(2*f-z*S.diff(f,z))+z*H*f)/B**3
 zero('source_A1_polynomial',A1.subs(sub,simultaneous=True)-P1)
 zero('source_A0_polynomial',A0.subs(sub,simultaneous=True)-P0)
 zero('source_A1_outer_linear',S.limit(r*A1,r,S.oo)-28*lam*m)
 zero('source_A0_outer_constant',S.limit(A0,r,S.oo)-56*lam**2*m)
 zero('potential_shift_outer_zero',S.limit(Fr*c,r,S.oo))
 zero('A1_zero_at_polynomial_horizon',P1.subs(x,-(f-x)))
 zero('A0_zero_at_polynomial_horizon',P0.subs(x,-(f-x)))
 zero('potential_shift_zero_at_polynomial_horizon',(f*cp/z**2).subs(x,-(f-x)))
 w=S.Function('w')(y);a=S.Function('a')(y);beta=S.Function('beta')(y)
 A1=S.Function('A1')(y);A0=S.Function('A0')(y)
 g=-A1*(S.diff(a,y)+S.diff(beta,y))+A0*beta
 Ka=S.diff(w*A1,y);Kb=Ka+w*A0
 zero('moment_adjoint_with_boundary_term',w*g-(Ka*a+Kb*beta-S.diff(w*A1*(a+beta),y)))
 V=S.Function('V')(y);U=S.Function('U')(y);gam=S.symbols('gamma')
 green=S.diff(S.diff(w,y)*U-w*S.diff(U,y),y)
 zero('projected_trace_Green_identity',green.subs(S.diff(w,y,2),(V+gam**2)*w)-w*(-S.diff(U,y,2)+(V+gam**2)*U))
 k,h0,Ac,As=S.symbols('k h0 Ac As',real=True)
 phase=Ac*S.cos(k*S.log(y))+As*S.sin(k*S.log(y))
 phaseT=-Ac*S.sin(k*S.log(y))+As*S.cos(k*S.log(y))
 ww=S.sqrt(y)*phase
 zero('outer_Ka_coefficient',S.diff(ww*h0*y,y)-h0*S.sqrt(y)*(S.Rational(3,2)*phase+k*phaseT))
 zero('outer_Kbeta_coefficient',S.diff(ww*h0*y,y)+2*h0*ww-h0*S.sqrt(y)*(S.Rational(7,2)*phase+k*phaseT))

class Dual:
 """Value and exact first variation, evaluated only at fixed rational jets."""
 def __init__(self,v=0,d=0):self.v=S.sympify(v);self.d=S.sympify(d)
 @staticmethod
 def coerce(x):return x if isinstance(x,Dual) else Dual(x)
 def __add__(self,x):
  x=self.coerce(x);return Dual(self.v+x.v,self.d+x.d)
 __radd__=__add__
 def __neg__(self):return Dual(-self.v,-self.d)
 def __sub__(self,x):return self+-self.coerce(x)
 def __rsub__(self,x):return self.coerce(x)+-self
 def __mul__(self,x):
  x=self.coerce(x);return Dual(self.v*x.v,self.d*x.v+self.v*x.d)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=self.coerce(x);return Dual(self.v/x.v,(self.d*x.v-self.v*x.d)/x.v**2)
 def __rtruediv__(self,x):return self.coerce(x)/self
 def __pow__(self,n):return Dual(self.v**n,n*self.v**(n-1)*self.d) if n else Dual(1)

def direct_curvature(F,a,b,t,r,point):
 th,ph=S.symbols('theta phi',real=True);coords=[t,r,th,ph]
 at={t:point[0],r:point[1],th:S.pi/2,ph:0}
 base=[-F,1/F,r*r,r*r*S.sin(th)**2]
 varied=[-b-2*F*a,-b/F**2,S.Integer(0),S.Integer(0)]
 def jet(i,*indices):
  v=base[i];dv=varied[i]
  for j in indices:v=S.diff(v,coords[j]);dv=S.diff(dv,coords[j])
  return Dual(v.subs(at),dv.subs(at))
 gg=[jet(i) for i in range(4)];inv=[1/g for g in gg]
 first={(j,i):jet(i,j) for j,i in product(range(4),repeat=2)}
 second={(k,j,i):jet(i,k,j) for k,j,i in product(range(4),repeat=3)}
 def dg(j,i,k):return first[j,i] if i==k else Dual()
 def ddg(l,j,i,k):return second[l,j,i] if i==k else Dual()
 def bracket(i,j,k):return dg(j,i,k)+dg(k,i,j)-dg(i,j,k)
 Gamma={(i,j,k):inv[i]*bracket(i,j,k)/2 for i,j,k in product(range(4),repeat=3)}
 def derivative(l,i,j,k):
  di=-inv[i]**2*first[l,i]
  return (di*bracket(i,j,k)+inv[i]*(ddg(l,j,i,k)+ddg(l,k,i,j)-ddg(l,i,j,k)))/2
 Riem={}
 for i,j,k,l in product(range(4),repeat=4):
  Riem[i,j,k,l]=derivative(k,i,l,j)-derivative(l,i,k,j)+sum(
   (Gamma[i,k,n]*Gamma[n,l,j]-Gamma[i,l,n]*Gamma[n,k,j] for n in range(4)),Dual())
 Ric={(j,l):sum((Riem[i,j,i,l] for i in range(4)),Dual()) for j,l in product(range(4),repeat=2)}
 scalar=sum((inv[i]*Ric[i,i] for i in range(4)),Dual())
 ric2=sum((inv[i]*inv[j]*Ric[i,j]**2 for i,j in product(range(4),repeat=2)),Dual())
 riem2=sum((gg[i]*inv[j]*inv[k]*inv[l]*Riem[i,j,k,l]**2 for i,j,k,l in product(range(4),repeat=4)),Dual())
 C2=riem2-2*ric2+scalar**2/3;J2=ric2-scalar**2/4
 Rbg=radial(F,r)[0]
 boxR=inv[1]*S.diff(Rbg,r,2).subs(at)-sum((inv[i]*Gamma[1,i,i] for i in range(4)),Dual())*S.diff(Rbg,r).subs(at)
 return scalar,C2,J2,boxR,at

def curvature_fixtures():
 t,r=S.symbols('t r',real=True)
 fixtures=[
  (1+2*r*r-6/r+3/r**2,t*r*r+t*t*r,t*t*r*r+t*r**3,(S.Rational(1,3),S.Integer(2))),
  (1+3*r*r/5-7/r+5/r**2-S.Rational(2,7)/r**3+S.Rational(1,7)/r**4,t*t/r+r*r,t*t*r+t/r**2+r**3,(S.Rational(2,5),S.Integer(3))),
  (1+2*r*r-6/r+3/r**2,t*r+r**3,t*t/r+t*r*r+r,(S.Rational(-1,4),S.Rational(5,2))),
  (1+2*r*r,t*t*r+r*r,t*t*r*r+t/r,(S.Rational(3,7),S.Integer(2)))]
 for ix,(F,a,b,pt) in enumerate(fixtures):
  scalar,C2,J2,boxR,at=direct_curvature(F,a,b,t,r,pt)
  R,W,j,psi,dC,dJ,KR,mix=predictions(F,a,b,t,r)
  ev=lambda x:S.simplify(x.subs(at));tag=str(ix)
  check('fixture_exterior_'+tag,ev(F)>0)
  zero('direct_4D_background_R_'+tag,scalar.v-ev(R))
  zero('direct_4D_background_C2_'+tag,C2.v-ev(W*W/3))
  zero('direct_4D_background_J2_'+tag,J2.v-ev(4*j*j))
  zero('direct_4D_delta_R_'+tag,scalar.d-ev(psi))
  zero('direct_4D_delta_C2_'+tag,C2.d-ev(dC))
  zero('direct_4D_delta_J2_'+tag,J2.d-ev(dJ))
  zero('direct_4D_delta_BoxR_'+tag,boxR.d-ev(KR))
  zero('direct_4D_complete_metric_coupling_'+tag,boxR.d+7*C2.d/12+11*J2.d/6-ev(mix))
  check('time_derivative_fixture_active_'+tag,ev(S.diff(b,t,2))!=0)
  NUMERICS.append({'fixture':'exact_4D_curvature','index':ix,'t':str(pt[0]),'r':str(pt[1]),
   'delta_R':str(scalar.d),'delta_C2':str(C2.d),'delta_J2':str(J2.d),'delta_BoxR':str(boxR.d)})

def moment_fixtures():
 r=S.symbols('r',positive=True);lam=S.Integer(2);ds=S.Rational(1,100);rh=S.Rational(11,10)
 ext=(r-1)**2/r**2*(1+lam*(r*r+2*r+3)+ds/r**2)
 F=S.factor(ext-rh*ext.subs(r,rh)/r)
 R,W,j,c,d,e=radial(F,r);A1=r*F*d;A0=r*F*(d*S.diff(F,r)+e*F)
 lo=S.Rational(13,10);hi=S.Rational(22,10);env=((r-lo)*(hi-r))**4
 a=env*(1+r);beta=env*(2-r)
 weights=[S.exp(-r),S.exp(-r)*(1+S.I*r/3),S.cos(2*r)/(1+r)]
 numeric=lambda ex:S.lambdify(r,ex,'mpmath')
 ff=numeric(F);aa1=numeric(A1);aa0=numeric(A0)
 for ix,w in enumerate(weights):
  # dy=-dr/F reverses endpoints; integrate from lo to hi with positive dr/F.
  gm=A1*F*(S.diff(a,r)+S.diff(beta,r))+A0*beta
  Ka=-F*S.diff(w*A1,r);Kb=Ka+w*A0
  direct=numeric(w*gm/F);adj=numeric((Ka*a+Kb*beta)/F)
  bounds=[mp.mpf('1.3'),mp.mpf('1.75'),mp.mpf('2.2')]
  val=mp.quad(direct,bounds);other=mp.quad(adj,bounds)
  error=abs(val-other)/(1+abs(val))
  check('adjoint_moment_quadrature_'+str(ix),error<mp.mpf(PROTOCOL['quadrature_relative_tolerance']),mp.nstr(error,15))
  check('test_moment_nonzero_'+str(ix),abs(val)>mp.mpf('1e-8'))
  NUMERICS.append({'fixture':'moment_identity_test_weight_not_horizon_solution','index':ix,
   'direct':mp.nstr(val,25),'adjoint':mp.nstr(other,25),'relative_error':mp.nstr(error,15)})
 z=S.symbols('z');fz=S.Poly(S.cancel(z*z*F.subs(r,1/z)),z)
 dz=S.Poly(S.cancel(d.subs(r,1/z)/z**4),z);ez=S.Poly(S.cancel(e.subs(r,1/z)/z**6),z)
 l1=lambda p:sum(abs(v) for v in p.all_coeffs())
 C1=l1(fz)*l1(dz);C0=l1(fz)*(l1(dz)*(2*l1(fz)+l1(S.Poly(z*S.diff(fz.as_expr(),z),z)))+l1(ez)*l1(fz))
 # Exact coefficient bounds are global; finite rational samples detect transcription errors.
 for ix,zz in enumerate([S.Rational(1,1000),S.Rational(1,10),S.Rational(1,2),S.Rational(9,10)]):
  rr=1/zz
  check('global_coefficient_bound_sample_'+str(ix),abs(A1.subs(r,rr))<=C1 and abs(A0.subs(r,rr))<=C0)
 check('fixture_horizon_simple',F.subs(r,rh)==0 and S.diff(F,r).subs(r,rh)>0)
 check('fixture_outer_source_nonzero',S.limit(A0,r,S.oo)>0)

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 mp.mp.dps=PROTOCOL['mpmath_precision'];ps=paths()
 hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 explicit metric source','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();curvature_fixtures();moment_fixtures()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 explicit metric source','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'exact_first_order_metric_source_derived':good,'global_H1_to_L2_source_bound_derived':good,
  'explicit_L2_moment_kernels_derived':good,'nontrivial_unconstrained_metric_moment_derived':good,
  'projected_trace_identity_retained':good,'all_local_contacts_retained':good,
  'physical_coupled_moment_evaluated':False,'physical_metric_boundary_selected':False,
  'metric_coupling_set_to_zero':False,'actual_metric_causal_inverse_proved':False,
  'actual_metric_mode_constructed':False,'actual_charged_stability_proved':False,
  'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,
  'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
