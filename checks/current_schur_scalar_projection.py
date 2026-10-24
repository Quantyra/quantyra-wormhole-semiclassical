"""Actual global scalar return and its source projection; shape row remains."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'checks'))
import current_schur_bolt_functional as prev
from flint import arb,arb_mat,ctx
import mpmath as mp
power=arb(0);y0=arb(0);moments={}
STEM='current-schur-scalar-projection';PREVIOUS='current-schur-bolt-functional'
OWN=['notes/'+STEM+'.md','checks/current_schur_scalar_projection.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=prev.encoded
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
prev.prior.check=check;prev.york.check=check
def digest(p):return sha256((root/p).read_bytes()).hexdigest()
def dec(x): return x.str(22)
class FP:
 def __init__(self,c=0):
  if isinstance(c,FP):self.c=dict(c.c)
  elif isinstance(c,dict):self.c={k:arb(v) for k,v in c.items() if not arb(v).is_zero()}
  else:self.c={} if arb(c).is_zero() else {(0,0,0):arb(c)}
 def __add__(self,b):
  b=FP(b);out=dict(self.c)
  for k,c in b.c.items():out[k]=out.get(k,arb(0))+c
  return FP(out)
 __radd__=__add__
 def __neg__(self):return FP({k:-c for k,c in self.c.items()})
 def __sub__(self,b):return self+-FP(b)
 def __rsub__(self,b):return FP(b)+-self
 def __mul__(self,b):
  b=FP(b);out={}
  for k,c in self.c.items():
   for l,d in b.c.items():
    key=tuple(x+y for x,y in zip(k,l));out[key]=out.get(key,arb(0))+c*d
  return FP(out)
 __rmul__=__mul__
 def __truediv__(self,b):return self*(1/arb(b))
 def __pow__(self,n):
  out=FP(1)
  for _ in range(n):out=out*self
  return out
 def diff(self):
  out=FP(0)
  for (k,j,t),c in self.c.items():
   if k:out+=FP({(k-1,j,t):k*c})
   if j or t:out+=FP({(k,j-1,t):-(j+t*power)*c})
  return out
 def integral(self):return sum((v*moment(*k) for k,v in self.c.items()),arb(0))
 def serial(self):return [{'z_power':k,'y_integer_power':j,'fractional_multiplicity':t,'coefficient':prev.encoded(c)} for (k,j,t),c in sorted(self.c.items())]
 def value(self,z):return sum((c*z**k*(1-z)**(j+t*power) for (k,j,t),c in self.c.items()),arb(0))
def mon(k=0,j=0,t=0):return FP({(k,j,t):1})
def moment(k,j,t):
 assert k>=0
 key=(k,j,t)
 if key not in moments:
  val=arb(0)
  for i in range(k+1):
   e=j+t*power+i+1
   v=-y0.log() if e.is_zero() else -(e*y0.log()).expm1()/e
   val+=((-1)**i)*comb(k,i)*v
  moments[key]=val
 return moments[key]
def norm(poly):
 mid=FP({k:v.mid() for k,v in poly.c.items()});sq=(mon(4)*mid**2).integral()
 assert sq>=0,dec(sq)
 err=sum((v.rad()*moment(4+2*k,2*j,2*t).sqrt() for (k,j,t),v in poly.c.items()),arb(0))
 return sq.sqrt()+err,sq,err
def bernstein(a):
 n=len(a)-1
 return [sum((a[i]*arb(comb(k,i))/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
def split(b):
 rows=[b]
 while len(rows[-1])>1:rows.append([(a+c)/2 for a,c in zip(rows[-1][:-1],rows[-1][1:])])
 return [r[0] for r in rows],[r[-1] for r in reversed(rows)]
def bern_product(a,b):
 m=len(a)-1;n=len(b)-1
 return [sum((a[i]*b[k-i]*comb(m,i)*comb(n,k-i)/comb(m+n,k)
  for i in range(max(0,k-n),min(m,k)+1)),arb(0)) for k in range(m+n+1)]
def taylor_lower(lo,hi,aa,bb):
 width=hi-lo;ymid=1-(hi+lo)/2;p=power-1;ww=[arb(0) for _ in range(7)];bc=arb(1)
 for j in range(7):
  if j:bc*= (p-j+1)/j
  coeff=ymid**p*bc*(-width/ymid)**j
  for i in range(j+1):ww[i]+=coeff*comb(j,i)*(-arb(1)/2)**(j-i)
 product=bern_product(bb,bernstein(ww));degree=len(product)-1
 raised=bern_product(aa,[arb(1)]*(degree-len(aa)+2))
 fall=arb(1)
 for j in range(7):fall*=p-j
 error=fall.abs_upper()/factorial(7)*(width/2)**7*(1-hi)**(p-7)*max(c.abs_upper() for c in bb)
 return min((a+b).lower() for a,b in zip(raised,product))-error
def monotonicity(v):
 n=max(k-4 for k,j,t in v.c)
 P=[v.c.get((i+4,0,0),arb(0)) for i in range(n+1)]
 T=[v.c.get((i+4,0,1),arb(0)) for i in range(n+1)]
 A=[(i+4)*c for i,c in enumerate(P)];B=[arb(0) for _ in range(n+2)]
 for i,c in enumerate(T):B[i]+=(i+4)*c;B[i+1]-=(i+4+power)*c
 cells=[];pending=[(arb(0),arb(1),bernstein(A),bernstein(B),0)]
 while pending:
  lo,hi,aa,bb,depth=pending.pop();al=min(x.lower() for x in aa);bl=min(x.lower() for x in bb)
  if bl>=0:lower=al+bl*(1-lo)**(power-1)
  elif hi<1:lower=al+bl*(1-hi)**(power-1)
  else:lower=arb(-1)
  method='direct_factor'
  if not lower>0 and hi<1:
   lower=taylor_lower(lo,hi,aa,bb);method='sixth_order_weight_Taylor_with_seventh_derivative_remainder'
  if lower>0:
   cells.append(dict(left=encoded(lo),right=encoded(hi),derivative_bracket_lower=encoded(lower),method=method))
  else:
   assert depth<24,(dec(lo),dec(hi),dec(al),dec(bl))
   am,ap=split(aa);bm,bp=split(bb);mid=(lo+hi)/2
   pending.extend([(mid,hi,ap,bp,depth+1),(lo,mid,am,bm,depth+1)])
 return cells
def convert(poly):return FP({(k,0,0):c for (k,l),c in poly.c.items() if l==0})
def absolute(poly):return FP({(k,0,0):c.abs_upper() for (k,l),c in poly.c.items() if l==0})
def execute(old):
 global power,y0
 p,h,zh,env=prev.york.actual_parameters();y0=h/(1+h);C=p['C'];B2=p['B2'];lam=p['g'];x=lam*B2
 NUMERICS.extend(prev.york.NUMERICS)
 power=((1+16*x/C).sqrt()-1).mid()/2
 beta=x/C;u=1/C;eps=p['eps'];v=1+h
 tau=v*h*h*(beta*(1+2/v+3/v**2)+u/v**2+eps/v**4)
 Q=mon(0,2)*(beta*(1+2*mon(1)+3*mon(2))+u*mon(2)+eps*mon(4))-tau*mon(3)
 R=-12*x+eps*C*(4*mon(5)-6*mon(6))
 tt=(1+3*x+eps*C)*mon(4)-eps*C*(5*mon(5)-arb(9)/2*mon(6))
 ee=1/(60*arb.pi());DR=96*arb.pi()*arb('1e22');D=DR-ee
 c4=6*D*C/B2;c1=6+(-2*DR+ee/6)*R/B2-22*ee/(3*B2)*tt
 mass=(2*R+48*arb(2000)**-10*B2)/C
 delta=lambda f:(Q*f.diff()).diff()
 op=lambda f:delta(f)-2*Q*mon(-1)*f.diff()
 adj=lambda f:delta(f)+(2*Q*mon(-1)*f).diff()
 strong=lambda f:c4*adj(op(f))-(c1*Q*mon(-2)*f.diff()).diff()+mass*mon(-4)*f
 weak=lambda a,b:(c4*op(a)*op(b)+c1*Q*mon(-2)*a.diff()*b.diff()+mass*mon(-4)*a*b).integral()
 gr=next(r for r in old['numerics'] if r['kind']=='actual_scalar_mixed_point_return' and r['name']=='TT_zero_J')
 g=FP({(c['power'],0,0):prev.endpoint.raw.decode(c['coefficient']) for c in gr['g_coefficients']})
 rhs=(x/C)*mon(-4)*g
 mstar=arb('3.99');mlo=mass.c[(0,0,0)].lower()-sum((v.abs_upper() for k,v in mass.c.items() if k!=(0,0,0)),arb(0))
 clo=c1.c[(0,0,0)].lower()-sum((v.abs_upper() for k,v in c1.c.items() if k!=(0,0,0)),arb(0))
 check('actual_scalar_mass_coercivity',mlo>mstar,dec(mlo))
 check('actual_scalar_gradient_coercivity',clo>arb('5.99'),dec(clo))
 check('actual_scalar_fourth_coefficient_positive',c4>0)
 check('fixed_positive_cold_exponent',power>arb('.45') and power<arb('.46'))
 # Independent polynomial weak assembly, using the already retained generator.
 LP=prev.local.Poly;lzp=prev.local.monomial
 Qp=LP({(i,0):env[s] for i,s in zip(prev.york.powers,prev.york.qs)})
 for k in [4,5]:
  known=sum(prev.local.scalar_covariant(p,Qp,lzp(k),lzp(4)).values(),LP(0)).integral(zh)
  check('independent_complete_weak_scalar_'+str(k),weak(mon(k),mon(4)).overlaps(known))
 for n in [11]:
  basis=[mon(i+4,0,t) for t in [0,1] for i in range(n)]
  # Choose one internally consistent point operator before assembling the
  # ill-conditioned trial matrix; restore full intervals for its residual.
  midpoly=lambda f:FP({k:c.mid() for k,c in f.c.items()})
  Qs=midpoly(Q);c1s=midpoly(c1);ms=midpoly(mass);rs=midpoly(rhs)
  ops=lambda f:(Qs*f.diff()).diff()-2*Qs*mon(-1)*f.diff()
  ws=lambda a,b:(c4.mid()*ops(a)*ops(b)+c1s*Qs*mon(-2)*a.diff()*b.diff()+ms*mon(-4)*a*b).integral()
  saved_y=y0;y0=y0.mid();moments.clear()
  mat=arb_mat([[ws(a,b) for b in basis] for a in basis])
  r=arb_mat([[(rs*a).integral()] for a in basis]);sol=mat.solve(r)
  trial=sum((f*sol[i,0].mid() for i,f in enumerate(basis)),FP(0));res=strong(trial)-rhs
  y0=saved_y;moments.clear()
  rn,sq,rad=norm(res)
  NUMERICS.append(dict(kind='scalar_trial_residual_before_acceptance',dimension=2*n,residual_weighted_L2_upper=encoded(rn),
   midpoint_residual_squared=encoded(sq),coefficient_uncertainty_L2_upper=encoded(rad),trial_coefficients=trial.serial()))
  check('whole_domain_scalar_residual_target',rn<arb('1e-9'),dec(rn))
  check('strong_weak_endpoint_agreement',((strong(trial)*mon(4)).integral()-weak(trial,mon(4))).contains(0))
  vtest=trial+mon(4)/2;cells=monotonicity(vtest)
  check('actual_trial_test_is_monotone',len(cells)>0,len(cells))
  # Exact known trace terms and complete upper allowances for unknown source contacts.
  du=p['eps']*C/x;eta=ee*lam;dl=D*lam;S=ee*arb(p['q'])**2/(lam*B2*B2)
  thermal_delta_over_x=tau*C/x
  fa=4+2/x+thermal_delta_over_x;fb=3+1/x+du
  ff=lzp(-2)+1/x-fa*lzp(1)+fb*lzp(2)-2*du*lzp(3)+du*lzp(4)
  Lap=lambda f:lzp(4)*(ff*f.derivative()).derivative()
  rr=4*lzp(5)-6*lzp(6);rnpoly=-12+du*rr
  tn=(3+1/x+du)*lzp(4)-du*(5*lzp(5)-arb(9)/2*lzp(6))
  # Cancel the exact outer constants algebraically, before interval majorants.
  cn=fa*lzp(3)-2*fb*lzp(4)+arb(20)/3*du*lzp(5)-5*du*lzp(6)
  check('correlated_outer_Weyl_identity',(cn-(-lzp(3)*ff.derivative()+rnpoly/6)).zero_enclosed())
  P=7*cn**2+arb(22)/3*tn**2-arb(11)/72*(rnpoly**2-144)
  Ps=2*(6+du)**2+arb(22)/3*(6+du)/x+arb(22)/3/x**2+22
  AA=-10*lzp(4)-20*lzp(5)+30*lzp(6)
  rem=eta*(P+Ps*rr/2)-dl*du*Lap(rr);tau0=S*AA+rem
  Kt=arb(11)/36*rnpoly-arb(7)/3*cn+arb(11)/3*tn
  K=dl*Lap(tau0)+eta*Kt*tau0
  Af=-14*eta*cn-dl*du*(lzp(2)*rr.derivative().derivative()+lzp(1)*rr.derivative())
  AQ=eta*(arb(44)/3*tn-14*cn)-dl*du*lzp(1)*rr.derivative()
  tmajor=-S*convert(AA)+absolute(rem)
  sr=arb('3.33e-10')*h;amax=arb('3.3e-10');H=arb('1.08e-9')
  # Keep each uncertain physical coefficient outside its exact trial moment.
  kval=(x/C)*sum((c*(mon(k-4)*vtest).integral() for (k,l),c in K.c.items()),arb(0))
  areal=(x/C)*(mon(-4)*(arb('2.11e-6')*absolute(Af)+arb('2.01e-7')*absolute(AQ))*vtest).integral()
  density=(x/C)*(amax*(mon(-4)*tmajor*vtest).integral()
   +sr*(mon(-2)*tmajor*vtest.diff()).integral()
   +2*H*(mon(-3,1)*tmajor*vtest.diff()).integral()
   +sr*zh**-2*tau0.value(zh).abs_upper()*vtest.value(zh))
  E=absolute(K)+arb('2.11e-6')*absolute(Af)+arb('2.01e-7')*absolute(AQ)+arb('8.7e-5')*absolute(tau0)
  E+=arb('2.493e-9')*absolute(lzp(1)*(1-lzp(1))*tau0.derivative())
  check('source_majorant_outer_integrability',all(k>=3 and j==t==0 for k,j,t in E.c))
  source_sq=(mon(-4)*E**2).integral();Cs=(x/C)*source_sq.sqrt()
  err=Cs*rn/mstar;allow=areal.upper()+density.upper()+err.upper()
  scalar_interval=kval+arb(0,allow)
  check('unevaluated_source_allowances_positive',areal>0 and density>0 and err>0)
  check('actual_scalar_inverse_pairing_error_target',err<arb('1e-40'),dec(err))
  threshold=arb('2.6035780e-40')
  NUMERICS.append(dict(kind='actual_global_scalar_inverse',dimension=2*n,h=encoded(h),zh=encoded(zh),y_lower=encoded(y0),
   fractional_power=encoded(power),c4=encoded(c4),gradient_lower=encoded(clo),mass_lower=encoded(mlo),mass_used=encoded(mstar),
   residual_weighted_L2_upper=encoded(rn),midpoint_residual_squared=encoded(sq),coefficient_uncertainty_L2_upper=encoded(rad),
   scalar_form_error_upper=encoded(rn/mstar.sqrt()),weighted_L2_error_upper=encoded(rn/mstar),
   trial_bolt_value=encoded(trial.value(zh)),trial_coefficients=trial.serial(),residual_coefficients=res.serial(),
   global_test_monotonicity_cells=cells,actual_full_scalar_operator_used=True,full_shape_source_evaluated=False))
  NUMERICS.append(dict(kind='actual_canonical_scalar_projection',action_normalization='R_nu(z^4/2+A^-1 B xi), in I/W units',
   known_areal_pairing=encoded(kval),unknown_areal_upper=encoded(areal),unknown_density_transport_endpoint_upper=encoded(density),
   source_mass_norm_upper=encoded(Cs),source_majorant=E.serial(),source_majorant_square_integral=encoded(source_sq),
   scalar_inverse_numerical_pairing_error=encoded(err),total_scalar_projection_enclosure=encoded(scalar_interval),
   necessary_canonical_shape_interval=encoded(scalar_interval+arb(0,threshold.upper())),
   source_uncertainty_to_reduced_threshold=encoded((areal+density)/threshold),
   scalar_numerical_error_to_reduced_threshold=encoded(err/threshold),
   full_signed_reduced_projection_computed=False,unknown_shape_may_cancel=True))
  print(json.dumps(dict(event='scalar_projection',scalar_enclosure=dec(scalar_interval),numerical_error=dec(err),
   unknown_contact_allowance=dec(areal+density),source_uncertainty_factor=dec((areal+density)/threshold))),flush=True)

def main():
 ctx.prec=512;mp.mp.dps=55
 review=json.loads((root/OWN[2]).read_bytes());old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
 before={p:digest(p) for p in paths};deps=prev.endpoint.dependencies();error=None
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],dependencies=deps,
  invocation=sys.argv,previous_goal_turn='progress')
 try:
  check('previous_actual_scalar_point_test_passed',old['status']=='passed')
  check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  check('pinned_interval_library',prev.prior.flint.__version__=='0.9.0')
  execute(old)
 except Exception:error=traceback.format_exc()
 after={p:digest(p) for p in paths};depafter=prev.endpoint.dependencies()
 status='passed' if error is None and before==after and deps==depafter else 'failed'
 run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
  checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
  actual_global_scalar_inverse_certified=status=='passed',canonical_scalar_source_projection_enclosed=status=='passed',
  signed_reduced_source_projection_executed=False,full_shape_source_evaluated=False,full_solution_feasibility_test_executed=False,
  nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,
  component_complete=False,goal_complete=False,formal_verification=False)
 out=root/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
  for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,
              dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
 print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
 if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
