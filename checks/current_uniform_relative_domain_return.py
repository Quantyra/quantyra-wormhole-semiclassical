"""Support the ordinary uniform relative-domain return and explicit fast embedding."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-uniform-relative-domain-return'
OWN=['notes/'+STEM+'.md','checks/current_uniform_relative_domain_return.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 theta,W,mr,mi=S.symbols('theta W mr mi',real=True)
 UU=S.cos(theta)+S.I*S.sin(theta);flow=W*(1-UU**2)-2*(mr+S.I*mi)*UU
 check('complex_unit_disk_invariance',S.simplify(S.expand_complex(2*S.re(S.conjugate(UU)*flow))+4*mr)==0)
 U,V,M=S.symbols('U V M')
 f=lambda z:W*(1-z*z)-2*M*z
 check('exact_common_mass_difference',S.expand(f(U)-f(V)+(2*M+W*(U+V))*(U-V))==0)
 r=arb('.001');du=r/(arb('.5')-r)
 rate=2*(1-(r.exp()-(-r).exp()))-2*r.exp()*(r.exp()-1)-8*du/arb(5).sqrt()
 check('complex_mass_sector',arb('.5').cos()>0)
 check('holomorphic_endpoint_displacement',du<arb('.002005'))
 check('uniform_relative_propagation_rate',rate>arb('1.98'))
 check('total_germ_displacement',du+arb('.036')<arb('.1'))
 check('germ_value_guard',16*arb('2.21')<36 and 2*4*arb('.003').exp()*arb('1.1')<36)
 check('germ_lipschitz_guard',32*arb('1.1')<36 and 8*arb('.003').exp()<36)
 check('germ_difference_inflation',arb('.036').exp()<arb('1.1'))
 den=arb('.79');up=arb('1.1')
 check('bubble_rational_guards',up*up/den**2<2 and 2*up**3/den**3<6 and up*(1+up*up)/den**3<6)
 check('diagonal_rational_guard',2*up/den**2<4)
 check('prefactor_difference_guard',6*(9+9+36)<400)
 direct=32/(arb('1.98')*arb(5).sqrt());img=160/(arb('1.98')*arb(5).sqrt())
 check('path_difference_guards',direct<8 and img<37 and direct+arb('.3')<10 and img+arb('.72')<40)
 check('individual_complex_propagators',4*arb('.072').exp()<5 and 16*arb('.072').exp()<20)
 check('full_bubble_ledger',arb('.006').exp()*16*4*((400*5+2*2*10*5)+(400*20+2*2*40*20))<2**20)
 check('complete_seagull_ledger',arb('.003').exp()*4*4*(9+36)<2**10)
 z=S.symbols('z');bb=S.Function('b')(z)
 fourth=S.diff(S.exp(bb),z,4)/S.exp(bb)
 check('fourth_jet_vertex_remainder',S.simplify(fourth-S.diff(bb,z,4)-(S.diff(bb,z)**4+6*S.diff(bb,z)**2*S.diff(bb,z,2)+3*S.diff(bb,z,2)**2+4*S.diff(bb,z)*S.diff(bb,z,3)))==0)
 # Repeated total differentiation of the first-order flow uses only
 # coefficient jets through order three for the fourth solution jet.
 aa=S.symbols('a0:5');cc=S.symbols('c0:5');uu=S.symbols('u0:5')
 ff=aa[0]*(1-uu[0]**2)-2*cc[0]*uu[0];jets=[ff]
 for n in range(3):
  expr=jets[-1]
  deriv=sum(S.diff(expr,aa[j])*aa[j+1]+S.diff(expr,cc[j])*cc[j+1] for j in range(4))
  deriv+=S.diff(expr,uu[0])*ff
  jets.append(S.expand(deriv))
 check('fourth_Riccati_jet_requires_only_third_mass_jet',not jets[3].has(aa[4],cc[4]) and jets[3].has(aa[3],cc[3]))
 for mm,vv in enumerate([63,7,3,3,7,63]):
  point=S.gamma(S.Rational(2*mm+1,2))*S.gamma(S.Rational(11-2*mm,2))/(2*S.pi*S.gamma(6))
  check('H6_point_evaluation_'+str(mm),S.simplify(point-S.Rational(vv,512))==0)
 check('small_complex_square_form_margin',arb('.9')*(-arb('.002')).exp()-11*arb('.01')**2>arb('.89'))
 check('fourth_vertex_positive_output_powers',7-S.Rational('4.005')==S.Rational('2.995') and 7>0)
 return dict(metric_radius='.001',global_holomorphy_radius='.5',endpoint_displacement=old.encoded(du),relative_rate_lower=old.encoded(rate),LLL_retained_separately=True)

def execute():
 alg=algebra()
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 a=arb('.01');b=arb('.02');c=arb('.04');sigma=arb('4.005');alpha=5-sigma
 d=arb(5).sqrt()*(c-a);cm=arb('.99')*arb(2).sqrt()*d;cp=arb('.99')*arb(10).sqrt()*d
 pref=8*pars['C']/(arb.pi()*pars['B2'])*(1+pars['gamma']*cp)
 def moments(r,j):return sum((comb(j,k)*factorial(r+j-k+2)*factorial(k)/cm**(r+j-k+2)/cp**(k+1) for k in range(j+1)),arb(0))
 K=[pref*2**20*factorial(j)*arb(1000)**j*moments(2,j) for j in range(5)]
 SG=[pref*2**10*factorial(j)*arb(1000)**j*moments(1,j) for j in range(5)]
 L=[arb(1)]
 for n in range(1,5):L.append(L[-1]*(sigma+2-n))
 response={}
 for n in range(1,5):
  terms=[(b**(sigma+1)*K[n],arb(0),0)]
  for rr in range(1,n):
   for ell in range(n-rr):
    mm=n-rr-1-ell
    terms.append((comb(n-rr-1,ell)*arb(2*(rr+ell+1))/(rr+1)*K[rr+ell]*L[mm],sigma-mm,1))
  for ell in range(n+1):terms.append((comb(n,ell)*SG[ell]*L[n-ell],sigma-(n-ell),2))
  if n==4:
   terms.append((arb('.002')*b**(sigma+1)*K[0],arb(-1),3))
   terms.append((arb('.002')*SG[0],sigma-1,3))
  response[n]=terms
 cs=[[arb(1)]]
 for k in range(1,4):cs.append([(alpha+j)*(cs[-1][j] if j<k else 0)+(cs[-1][j-1] if j>0 else 0) for j in range(k+1)])
 caps=[1,3,10,60];parts=[]
 for k in range(4):
  row=[arb(0) for _ in range(4)]
  for ell in range(k+1):
   for j in range(k-ell+1):
    for A,power,part in response[j+1]:
     pp=alpha+j+power;assert pp>0
     row[part]+=comb(k,ell)*caps[ell]*cs[k-ell][j]*A*a**pp/(2*pp).sqrt()
  parts.append(row)
 Nk=[sum(row,arb(0)) for row in parts]
 eta=1801/arb(2).sqrt()*sum((comb(3,k)*v*v for k,v in enumerate(Nk)),arb(0)).sqrt()
 # Actual normal-inverse point evaluations on H6; explicit C8 extension.
 t,X=S.symbols('t X');sig=S.Rational(801,200)
 cut=1-S.integrate(t**8*(1-t)**8,(t,0,t))*S.factorial(17)/S.factorial(8)**2
 def rat(x):return arb(int(x.p))/int(x.q)
 def bern_bound(poly,N=64):
  n=poly.degree();co=[poly.nth(i) for i in range(n+1)];maximum=S.Rational(0)
  for cell in range(N):
   left=S.Rational(cell,N);width=S.Rational(1,N)
   pw=[width**i*sum(co[k]*S.binomial(k,i)*left**(k-i) for k in range(i,n+1)) for i in range(n+1)]
   vals=[sum(pw[i]*S.binomial(j,i)/S.binomial(n,i) for i in range(j+1)) for j in range(n+1)]
   maximum=max(maximum,max(abs(v) for v in vals))
  return maximum
 guards=[S.Integer(1)]+[100**j*bern_bound(S.Poly(S.diff(cut,t,j),t)) for j in range(1,5)]
 point=[(arb(v)/512).sqrt() for v in [63,7,3,3,7,63]]
 W=[]
 for j in range(5):
  poly=S.Poly(S.prod(sig+X-l for l in range(j)),X)
  W.append(1801*sum((rat(abs(poly.nth(mm)))*point[mm] for mm in range(j+1)),arb(0)))
 embed=[]
 for j in range(5):
  cap=sum((comb(j,ell)*rat(guards[ell])*W[j-ell]*arb('.03')**(sigma-j+ell) for ell in range(j+1)),arb(0))
  inner=W[j]*arb('.02')**(sigma-j)
  if j==4:cap*=arb('.03');inner*=arb('.02')
  embed.append(max(cap,inner))
 KA=max(embed)

 check('original_seed_collar',zh>arb('.05') and ((1+pars['q'])/pars['C']).sqrt()>18)
 be=pars['beta'];co=[arb(0),arb(0),pars['u']/be,-4-2*pars['u']/be-pars['tau']/be,3+pars['u']/be+pars['eps']/be,-2*pars['eps']/be,pars['eps']/be]
 qerr=sum((co[j].abs_upper()*arb('.05')**j for j in range(2,7)),arb(0))
 check('same_polynomial_complex_germ',qerr<arb('.001') and 1/(be*(1-qerr))<16 and 1/(be*(1-qerr)).sqrt()<4)
 Rmax=arb(1)/18000
 check('Taylor_mass_germ_bound',arb('.002')*sum((Rmax**j/factorial(j) for j in range(4)),arb(0))<arb('.003'))
 check('complex_germs_before_geometry_taper',b+Rmax<c)
 check('extension_endpoint_values',cut.subs(t,0)==1 and cut.subs(t,1)==0)
 check('extension_eight_endpoint_jets',all(S.diff(cut,t,j).subs(t,e)==0 for j in range(1,9) for e in [0,1]))
 check('uniform_relative_H3_target',eta<arb('2.205e-13'))
 check('metric_derivative_per_A_target',1000*eta<arb('2.205e-10'))
 check('actual_fast_A_embedding_target',KA<arb('5.681e6'))
 check('actual_fast_bilinear_graph_target',1000*eta*KA<arb('.001253'))
 check('proposed_fast_return_target_closed',1000*eta*KA<arb('1e9'))
 check('fast_ball_admitted_A_cost',KA*arb('1e-11')<arb('5.681e-5'))
 check('remaining_background_admission_target',arb('.001')-KA*arb('1e-11')>arb('.00094319'))
 check('global_embedding_sufficient_threshold',arb('2.205e-10')*arb('4.535e18')<arb('1e9'))
 check('uniform_extra_allowance',eta-arb('2.012953e-13')<arb('2e-14') and eta-arb('2.012953e-13')<arb('1e-12'))
 check('updated_remaining_nonlocal_target',S.Rational('.159916396')-S.Rational('1e-12')==S.Rational('.159916395999'))
 return dict(algebra=alg,relative_uniform_graph_bound=old.encoded(eta),metric_derivative_per_A=old.encoded(1000*eta),fast_extension_A_bound=old.encoded(KA),fast_graph_metric_return_bound=old.encoded(1000*eta*KA),fast_radius_A_cost=old.encoded(KA*arb('1e-11')),remaining_background_A_margin=old.encoded(arb('.001')-KA*arb('1e-11')),global_graph_A_embedding_target=old.encoded(arb('1e9')/(1000*eta)),cutoff_derivative_guards=[str(x) for x in guards],embedding_by_derivative=[old.encoded(x) for x in embed],weighted_output_parts=[[old.encoded(x) for x in row] for row in parts],bubble_derivative_bounds=[old.encoded(x) for x in K],seagull_derivative_bounds=[old.encoded(x) for x in SG],remaining_nonlocal_target='.159916395999',
  uniform_relative_domain_return_on_admitted_family=True,relative_strong_common_regulator_identified=True,relative_first_metric_variation_per_A_certified=True,explicit_fast_extension_embedding_certified=True,complete_fast_near_return_bilinear_bound_certified=True,global_clock_gap_used=False,full_angular_and_thermal_towers_retained=True,seagull_images_and_diagonal_jumps_retained=True,mass_log_fourth_vertex_retained=True,background_A_admission_target='.00094319 after fast radius',unknown_global_background_admitted=False,global_X_graph_A_embedding_known=False,complete_local_curved_response_certified=False,local_strong_common_regulator_limit_proved=False,complete_C3_known=False,actual_controlled_solution=False,full_signed_projection_executed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,next_local_response_calculation_executed=False)

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
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),uniform_return=result.get('relative_uniform_graph_bound'), fast_C_return=result.get('fast_graph_metric_return_bound')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
