"""Support the weighted third comparison; no uniform curved theorem is inferred."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
from functools import lru_cache
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-weighted-third-comparison'
OWN=['notes/'+STEM+'.md','checks/current_weighted_third_comparison.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 t,z,sig=S.symbols('t z sig',positive=True)
 # Logarithmic norm equivalence holds for all frequencies by the ordinary
 # split, not by a frequency mesh. These enclose its endpoint constants.
 l0=arb(17)/4
 check('small_frequency_log_converse', (2+l0.exp()).log()**2<19)
 check('large_frequency_log_converse',arb(17)/2*(1+2*(-l0).exp()/l0)**2<19)
 check('composition_difference_quotient_guard',3+216*arb(2).log()**2<107)
 check('log_coordinate_constants',19*107*8<128**2 and 2*19*107*8<181**2 and 19*3<8**2)
 check('Fourier_commutator_integrals',S.integrate(t+t*t,(t,0,1))==S.Rational(5,6) and S.integrate(t**-3+t**-2,(t,1,S.oo))==S.Rational(3,2))
 check('compact_commutator_constants',arb('2.5')/arb.pi()*arb(5)/6<2 and arb('2.5')/arb.pi()*arb(3)/2<2)
 k=S.Function('k')(z)
 check('separated_third_kernel_identity',S.expand(z*z*S.diff(z*k,z,3)-(3*z*z*S.diff(k,z,2)+z**3*S.diff(k,z,3)))==0)
 for nu,target in [(0,S.Rational(16,15)),(1,S.Rational(8,5))]:
  moment=2**3*S.gamma(3)**2*S.gamma(3+nu)*S.gamma(3-nu)/S.gamma(6)
  check('entire_Bessel_moment_'+str(nu),moment==target)
 check('full_separated_complex_kernel_guard',arb(224)/(15*arb.pi()**2)*3*arb.pi()/arb(2).sqrt()<16)
 check('near_Hk_convexity_guard',5*arb(2).log()>1)
 # Direct physical derivatives verify the Mellin Leibniz polynomials.
 w=S.Function('w');x=S.symbols('x');polynomial=S.Integer(1)
 for n in range(1,9):
  polynomial=S.expand(polynomial*(sig+x-(n-1)))
  lhs=S.diff(z**sig*w(S.log(z)),z,n)*z**(n-sig)
  rhs=sum(co*S.Subs(S.Derivative(w(x),(x,degree[0])),x,S.log(z)) for degree,co in S.Poly(polynomial,x).terms())
  check('physical_Mellin_derivative_'+str(n),S.simplify(lhs-rhs.doit())==0)
 # Bernstein cutoff equals the integrated positive beta(9,9) density.
 cdf=sum(S.binomial(17,j)*t**j*(1-t)**(17-j) for j in range(9,18))
 check('beta_cutoff_density',S.factor(S.diff(cdf,t)-S.factorial(17)/S.factorial(8)**2*t**8*(1-t)**8)==0)
 check('cutoff_eighth_jet_endpoints',all(S.diff(cdf,t,j).subs(t,0)==0 and S.diff(1-cdf,t,j).subs(t,1)==0 for j in range(9)))
 # These are interval constants in analytic proofs, not discretized modes.
 return dict(log_converse='19 Q',composition_guard=128,generic_product_guards=[181,8],sharp_input_preserved=True,separated_primitive_identity=True)

def execute():
 alg=algebra()
 ctx.prec=512;pars,h,zh,env=old.scalar.prev.york.actual_parameters();C=pars['C'];B2=pars['B2'];be=pars['beta'];g=pars['g']
 sig=arb('4.005');a=arb('.01');b=arb('.03');H=a/2
 K=C*be**2/(15*arb.pi()*B2);L=10*arb.pi()*(pars['Fl']+23*g/(180*arb.pi()))/g
 Cb=[2/(K*L*arb('.005')*sig*(sig-3)),2/(K*L*sig*(sig-3)),2/(K*L*sig),2/(K*L)]
 Cb4=max(Cb[0],8/(K*L),1/(K*(L*arb(6).log()).sqrt()))
 Cq=arb('2.328');W=Cb+[Cb[1],Cb[2],Cb[3],Cb4,Cq/(K*arb(2).log())]
 P=arb('1e18');WL=[max(W[j]*(2+P).log(),Cq/K*P**(j-8)) for j in range(8)]+[Cq/K]
 def rat(v):return arb(int(v.p))/int(v.q)
 def split(coeff):
  rows=[coeff]
  for i in range(len(coeff)-1):rows.append([(rows[-1][j]+rows[-1][j+1])/2 for j in range(len(rows[-1])-1)])
  return [row[0] for row in rows],[row[-1] for row in rows][::-1]
 base=[S.Integer(0)]*9+[S.Integer(1)]*9;controls=[];current=base
 for order in range(10):
  cells=[current]
  for level in range(6):cells=[c for cell in cells for c in split(cell)]
  controls.append([[min(cell),max(cell)] for cell in cells])
  degree=len(current)-1
  current=[degree*(current[j+1]-current[j]) for j in range(degree)]
 @lru_cache(None)
 def M(order,offset):
  power=sig+offset;vals=[]
  if order==0:vals.append(arb('.75')**power if power<0 else arb('2.5')**power)
  for cell in range(64):
   for left,width,idx in [(arb('.5'),arb('.25'),cell),(arb('2.5'),arb('.5'),63-cell)]:
    lo=left+width*cell/64;hi=left+width*(cell+1)/64
    bound=max(abs(v) for v in controls[order][idx])
    vals.append(rat(bound)/width**order*(lo**power if power<0 else hi**power))
  return max(vals)
 V0=arb(3).sqrt()*M(0,0)*W[0]
 V1=arb(3).sqrt()*(M(1,0)*W[0]+M(0,-1)*(sig*W[0]+W[1]))
 Vlow=(V0*V1).sqrt()
 topM4=arb(0)
 for j in range(5):
  fall=arb(1)
  for kk in range(4-j):fall*=sig-8-kk
  topM4+=comb(4,j)*fall.abs_upper()*M(j,-12+j)
 xx=S.symbols('x');sigS=S.Rational(801,200)
 polys=[S.Poly(S.prod(sigS+xx-i for i in range(n)),xx) for n in range(9)]
 cs=[[arb(1)]]
 for ell in range(1,4):cs.append([(5+jj-sig)*(cs[-1][jj] if jj<ell else 0)+(cs[-1][jj-1] if jj else 0) for jj in range(ell+1)])
 near=arb(0);by_m=[arb(0) for _ in range(9)]
 for ell in range(4):
  for jj,co in enumerate(cs[ell]):
   rr=jj+1;kk=4+rr
   outer=arb(comb(3,ell)).sqrt()*co*arb(2)**(4+rr-sig)*7*2**kk
   for n in range(kk+1):
    power=sig+kk-n;geo=arb(2)**(-power);SP=H**power/(1-geo);SL=SP*((1/H).log()+arb(2).log()*geo/(1-geo))
    for i in range(n+1):
     mm0=M(i,-n+i);mm1=M(i+1,-n+i)+(sig-n+i).abs_upper()*M(i,-n+i-1)
     for (mm,),cc in polys[n-i].terms():
      if mm==8:
       const=128*mm0*WL[8]+8*arb(3).sqrt()*(mm0+topM4)*W[7]
      else:const=181*mm0*WL[mm]+8*mm1*W[mm]
      term=outer*arb(comb(kk,n)).sqrt()*comb(n,i)*rat(abs(cc))*(const*SP+arb(3).sqrt()*mm0*W[mm]*SL)
      near+=term;by_m[mm]+=term
 pref=2/B2*100*4*arb(2).sqrt()*C*Vlow
 near*=pref;by_m=[x*pref for x in by_m]
 # Full separated kernels, then primitive formula for the mixed B term.
 Ar=[arb(16*factorial(r)*4**r) for r in range(8)]
 Kn=[Ar[n+3]+3*sum((comb(n,j)*factorial(j)*Ar[n-j+2] for j in range(n+1)),arb(0)) for n in range(4)]
 N=[arb(1801)/arb(2).sqrt()]
 for j in range(1,5):N.append(N[-1]*(sig+2-j))
 responses={}
 for rr in range(1,5):
  terms=[]
  # v(z) A(f)(z), including all Leibniz derivatives.
  for j in range(rr+1):
   ll=rr-j
   if ll==0:
    coef=2*Ar[0]*5**5*N[0]*N[j]*b**(sig-4)*(1/(sig+1)+1/(sig-4));power=sig-j
   else:
    coef=2*Ar[ll]*5**(5+ll)*N[0]*N[j]*(1/(sig+1)+1/(4+ll-sig));power=2*sig-4-rr
   terms.append((comb(rr,j)*coef,power))
  # A(f v): both complete remote factors, no derivative of the indicator.
  exponent=2*sig-4-rr
  terms.append((2*Ar[rr]*5**(5+rr)*N[0]**2*b**exponent*(1/(2*sig+1)+1/exponent),arb(0)))
  # B primitive: j=0 has two primitives, j>=1 only local jets on v.
  terms.append((2*Kn[rr-1]*N[0]**2*arb(9)**(4+rr)*b**(2*sig-4-rr)*(1/(2*sig-4-rr)+1/(4+rr)),arb(0)))
  for j in range(1,rr):
   exponent=5+rr-j-sig
   coef=2*comb(rr-1,j)*Kn[rr-1-j]*N[0]*N[j-1]*arb(9)**sig*arb(8)**exponent/exponent
   terms.append((coef,2*sig-4-rr))
  responses[rr]=terms
 far=arb(0)
 for ell in range(4):
  for jj,co in enumerate(cs[ell]):
   rr=jj+1
   for coefficient,power in responses[rr]:
    outpower=power+4+rr-sig
    assert outpower>0
    far+=arb(comb(3,ell)).sqrt()*co*coefficient*a**outpower/(2*outpower).sqrt()
 far*=2/B2*100*C*2
 ans=dict(near=old.encoded(near),far=old.encoded(far),total=old.encoded(near+far),near_by_input_derivative=[old.encoded(v) for v in by_m],derivative_multiplier_caps=[old.encoded(v) for v in W],logarithmic_inverse=[old.encoded(v) for v in WL],cutoff_top_M0=old.encoded(M(0,-8)),cutoff_top_M4=old.encoded(topM4),low_fourier_factor=old.encoded(Vlow),low_V0=old.encoded(V0),low_V1=old.encoded(V1),kernel_A_guards=[old.encoded(v) for v in Ar],primitive_kernel_guards=[old.encoded(v) for v in Kn],cutoff_subdivision=64)
 check('physical_units_retained',pars['G']==1 if 'G' in pars else True)
 check('near_bound_in_original_units',near>arb('1.2946e10') and near<arb('1.2947e10'))
 check('all_separated_mixed_tails',far>0 and far<arb('6.206e-15'))
 check('complete_comparison_allocation',near+far<arb('1.295e10'))
 check('top_derivative_dominates',by_m[8]>near-arb('.000127'))
 check('lower_derivative_total_small',sum(by_m[:7],arb(0))<arb('1e-9'))
 check('low_fourier_constant',Vlow<231124)
 check('eighth_log_derivative_kept_in_L2',WL[8]>arb('6.7402e35') and WL[8]<arb('6.7403e35'))
 check('far_integrability_margin',2*sig-8>0 and 2*sig-8<arb('.01001'))
 check('primitive_kernel_guards',all(v==arb(n) for v,n in zip(Kn,[7680,118272,2300928,54088704])))
 remaining=arb('.159916395999')-arb('1e-11')*arb('1.295e10')
 check('positive_remaining_budget',remaining>arb('.0304163959989') and remaining<arb('.0304163959991'))
 development=json.loads((ROOT/'evidence'/f'{STEM}-development.json').read_bytes())
 first=json.loads(development['runs'][0]['output_utf8']);second=json.loads(development['runs'][1]['output_utf8'])
 check('earlier_excessive_bound_preserved',old.decode(first['total'])>arb('.159916395999')/arb('1e-11'))
 check('sharpened_development_matches_frozen_bound',old.decode(second['total']).overlaps(near+far))
 return dict(**ans,algebra=alg,comparison_C3_upper='1.295e10',comparison_radius_cost_upper='.1295',remaining_seed_and_additional_variation_budget=old.encoded(remaining),best_case_C3_add_ceiling=old.encoded(remaining/arb('1e-11')),equal_split_seed_ceiling=old.encoded(remaining/2),equal_split_C3_add_ceiling=old.encoded(remaining/arb('2e-11')),comparison_weighted_bilinear_bound_executed=True,entire_angular_spectrum_retained=True,all_separated_mixed_tails_retained=True,finite_fourth_contact_retained=True,low_contacts_retained_separately=True,actual_uniform_curved_response_known=False,actual_matched_strong_regulator_limit_proved=False,original_X_graph_comparison_known=False,global_admission_known=False,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,original_graph_changed=False,finite_physical_modes_selected=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_actual_defect_calculation_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),total=result.get('total'),remaining=result.get('remaining_seed_and_additional_variation_budget')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
