"""Actual finite-geometry prescribed magnetic reference response, not h1."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
import sympy as sp
from flint import arb,ctx
STEM='current-magnetic-reference-reconstruction'
OWN=['notes/'+STEM+'.md','checks/current_magnetic_reference_reconstruction.py','evidence/'+STEM+'-review.json']
ROWS=[];enc=old.encoded;dec=old.decode
def check(name,ok,detail=''):
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def hull(lo,hi):return arb((lo+hi)/2,(hi-lo)/2)
def g(z):
 if not z>arb('.5'):return old.g(z)
 y=1-z;value=-4*y*y;yp=y*y
 for j in range(3,49):
  yp*=y;value-=2*yp/j
 return value+arb(0,(2*yp*y/(49*(1-y))).abs_upper())
def A(z):
 if not z>arb('.5'):return (3-2*z.log())*z-2*z*z+z**3-2
 y=1-z;value=-arb(4)/3*y**3;yp=y**3
 for j in range(3,49):
  yp*=y;value-=2*yp/(j*(j+1))
 return value+arb(0,(2*yp*y/(49*50*(1-y))).abs_upper())
def Q(p,z):return (1-z)**2*(p['beta']*(1+2*z+3*z*z)+p['u']*z*z+p['eps']*z**4)-p['tau']*z**3
def algebra():
 t,z,k,C,ga,vh,ss=sp.symbols('t z k C gamma v sigma',positive=True)
 gg=4/t-3/t**2-1-2*sp.log(t);aa=(3-2*sp.log(z))*z-2*z*z+z**3-2
 check('primitive_identity',sp.simplify(sp.diff(aa,z)+gg.subs(t,1/z))==0)
 check('source_double_product_zero',gg.subs(t,1)==0 and sp.diff(gg,t).subs(t,1)==0)
 check('signed_source_derivative',sp.simplify(sp.diff(gg,t)+2*(t-1)*(t+3)/t**3)==0)
 f=(-2*C*ga*vh*ss+k*C*(aa.subs(z,1/t)-aa.subs(z,1/vh)))/t
 check('first_reference_constraint',sp.simplify(sp.diff(t*f,t)-k*C*gg/t**2)==0)
 check('moving_bolt_constraint',sp.simplify(f.subs(t,vh)+2*C*ga*ss)==0)
 ftt=sp.symbols('Ftt')
 sig=-k*C*gg.subs(t,vh)/(vh**2*(vh*ftt+2*C*ga))
 check('fixed_temperature_constraint',sp.simplify((sp.diff(f,t).subs(t,vh)+ftt*ss).subs(ss,sig))==0)
 ww,ff=sp.Function('w')(t),sp.Function('F')(t)
 # Differentiating -F times the outward Green integral gives this ODE.
 check('fixed_chart_displacement_ODE',sp.simplify(sp.diff(-ff*sp.Integral(f/ff**2,(t,t,sp.oo)),t)-(-sp.diff(ff,t)*sp.Integral(f/ff**2,(t,t,sp.oo))+f/ff))==0)

def execute():
 algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 C,B2,x,ga=p['C'],p['B2'],p['g']*p['B2'],p['gamma'];vh=1+h
 k=arb(p['q'])**2/(12*arb.pi()*C*B2);Fl=p['Fl'];zo=arb(2)**-20;yc=arb(2)**-49;zc=1-yc
 z=sp.symbols('z');qexpr=(1-z)**2*(sp.Symbol('b')*(1+2*z+3*z*z)+sp.Symbol('u')*z*z+sp.Symbol('e')*z**4)-sp.Symbol('ta')*z**3
 # Exact coefficient evaluation; use Q(zh)=0 and Q'(zh)=-2 gamma.
 q2=2*p['u']+6*(-4*p['beta']-2*p['u']-p['tau'])*zh+12*(3*p['beta']+p['u']+p['eps'])*zh**2-40*p['eps']*zh**3+30*p['eps']*zh**4
 sigbar=-g(zh)*zh**2/(6*ga+zh*q2)
 check('actual_positive_bolt_shift',sigbar>0 and sigbar<3*h*h)
 Hbar=(8+15*h*h)*C/(arb('.99')*x)
 check('global_reference_shape_bound',Hbar<49)
 check('inner_regular_H_bound',(4+arb('7.5')*h*h)*h/ga<arb('4.1'))
 tail_log=(1+2*C*ga/(arb('.95')*x*h)).log()
 check('inner_outward_Green_tail',tail_log<arb('2.8'))
 Jbar=(arb('4.1')+arb('2.8')*49)/(2*C*ga)
 b_inner=arb(49)/2+3*h*h+((2*C*ga+arb('12.1')*x*h)/2+(2*C*ga*h+arb('6.05')*x*h*h))*Jbar
 check('inner_shape_amplitude',b_inner<170)
 check('outer_shape_amplitude',(arb('.5')+10)*49+3*h*h<515)
 check('outer_Q_guard',Q(p,hull(arb(0),zo))>p['beta']/2)
 seeddata=json.loads((ROOT/'evidence/current-schur-scalar-projection.json').read_bytes())['runs'][-1]
 row=next(r for r in seeddata['numerics'] if r['kind']=='actual_global_scalar_inverse')
 old.scalar.power=dec(row['fractional_power']);old.scalar.y0=h/(1+h)
 trial=old.scalar.FP({(r['z_power'],r['y_integer_power'],r['fractional_multiplicity']):dec(r['coefficient']) for r in row['trial_coefficients']})+old.scalar.mon(4)/2
 power=old.scalar.power;n=max(i for i,j,t in trial.c)
 P=[trial.c.get((i,0,0),arb(0)) for i in range(n+1)];T=[trial.c.get((i,0,1),arb(0)) for i in range(n+1)]
 AP=[arb(0) for _ in P];BT=[arb(0) for _ in P]
 for i in range(n+1):
  if i:AP[i-1]+=i*P[i];BT[i-1]+=i*T[i]
  AP[i]-=i*P[i];BT[i]-=(i+power)*T[i]
 C1=sum((a.abs_upper() for a in AP+BT),arb(0))
 initial=[(arb(2)**-j,arb(2)**(-j+1)) for j in range(20,1,-1)]
 initial +=[(1-arb(2)**(-j+1),1-arb(2)**-j) for j in range(2,50)]
 shift=A(zh)+2*ga*vh*sigbar
 check('outer_numerator_guard',2+shift.abs_upper()<3)
 Slo=arb(0);Shi=3*zo**4/p['beta']**2;integ=arb(0);cells=[];sub=128
 for band,(left,right) in enumerate(initial):
  for j in range(sub):
   lo=left+(right-left)*j/sub;hi=left+(right-left)*(j+1)/sub;zz=hull(lo,hi);width=hi-lo
   qq=Q(p,zz);check_positive=qq.lower()>0
   if not check_positive:raise AssertionError('positive_full_Q_'+str(band)+'_'+str(j))
   aa=hull(A(lo).lower(),A(hi).upper());num=-aa+shift
   if not num.lower()>0:raise AssertionError('positive_full_integrand_'+str(band)+'_'+str(j))
   ii=zz**3*num/qq**2
   Srange=hull(Slo.lower(),(Shi+width*ii.upper()).upper())
   wbar=qq/zz**2*Srange
   bern=[old.bern_on_cell(c,lo,hi) for c in [P,T,AP,BT]]
   vv=old.family_bound(bern[0],bern[1],lo,hi,power)
   yvp=old.family_bound(bern[2],bern[3],lo,hi,power)
   gg=hull(g(lo).lower(),g(hi).upper())
   kw=-2*zz**2*(1+3*zz)*yvp-8*gg*zz**5
   val=wbar*kw*width;integ+=val
   Slo+=width*ii.lower();Shi+=width*ii.upper()
   cells.append(dict(lo=enc(lo),hi=enc(hi),inner_integrand=enc(ii),w_over_k=enc(wbar),signed_integral=enc(val)))
  if band%12==0:print(json.dumps(dict(event='actual_nested_integral_band',band=band,cells=len(cells),integral=str(integ))),flush=True)
 check('whole_bulk_coverage',len(cells)==len(initial)*sub and initial[0][0]==zo and initial[-1][1]==zc)
 wend=Q(p,zc)/zc**2*hull(Slo.lower(),Shi.upper())
 gp=-2*yc*(1+3*zc);kb=2*gp*trial.value(zc)-4*g(zc)*zc**4
 endpoint=-zc**2*kb*wend/2
 logm=zo**5*(arb(1)/25-zo.log()/5)
 outer=98*((8*C1+64)*zo**5/5+16*logm)+3*h*h*((8*C1+64)*zo**6/6+16*zo**6*(arb(1)/36-zo.log()/6))
 n_cold=3*h*h+98*yc/zc
 kn=arb('8.08')+16*yc;kbmax=arb('16.16')*yc+16*yc*yc/zc
 cold=yc*(kn*n_cold+kbmax*600)
 scalarE=4392/arb(5).sqrt()*arb('3.6891662e-12')
 value=k*k/Fl*(integ+endpoint);error=k*k/Fl*(outer+cold+scalarE)
 enclosure=value+arb(0,error.upper())
 olddata=json.loads((ROOT/'evidence/current-product-matched-source-projection.json').read_bytes())['runs'][-1]
 r=next(a for a in olddata['numerics'] if a['kind']=='validated_magnetic_weak_integrals')
 iw=dec(r['integral_enclosures'][1]);raw=dec(r['raw_absolute_kernel_integral_bounds'][1])
 direct=k*k/Fl*(98*iw.upper()+3*h*h*raw+endpoint.abs_upper()+cold+scalarE)
 ceiling=arb('1.06654e-54')
 check('independent_unsigned_envelope_target',direct<arb('3e-57'))
 check('signed_complete_component_target',enclosure.abs_upper()<arb('3e-57'))
 check('component_enclosures_consistent',enclosure.overlaps(arb(0,direct.upper())))
 check('physical_bolt_component_target',k*sigbar/h<arb('3e-62'))
 remainder_target=ceiling-enclosure.abs_upper()
 return dict(parameters={a:enc(b) for a,b in dict(k_M=k,h=h,zh=zh,gamma=ga,C=C,B2=B2,F_lambda=Fl).items()},
  sigma_over_k=enc(sigbar),physical_bolt_component=enc(k*sigbar/h),global_nu_bound=enc(99*k),global_b_bound=enc(600*k),
  first_reference_shape_bound_over_k=enc(Hbar),inner_shape_bound_over_k=enc(b_inner),
  signed_bulk_coefficient=enc(integ),retained_endpoint_coefficient=enc(endpoint),
  normalized_projection_value=enc(value),outer_error=enc(k*k/Fl*outer),cold_error=enc(k*k/Fl*cold),
  scalar_trial_error=enc(k*k/Fl*scalarE),complete_signed_component_interval=enc(enclosure),
  independent_unsigned_component_upper=enc(direct),old_localized_necessary_ceiling_lower=enc(ceiling),
  sufficient_cross_projection_target=enc(remainder_target),ceiling_to_component_allowance=enc(ceiling/enclosure.abs_upper()),
  cells=cells,actual_component_reconstructed=True,full_first_correction_reconstructed=False,
  full_quantum_source_reconstructed=False,complete_cross_projection_known=False,actual_controlled_solution=False,
  nonlinear_admission=False,all_sector_stability=False)

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
 result=run.get('result',{});print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),summary={k:v for k,v in result.items() if k!='cells'}),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
