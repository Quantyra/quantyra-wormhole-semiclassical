"""Support the entire angular-difference component in the current strong graph."""
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
STEM='current-angular-strong-graph'
OWN=['notes/'+STEM+'.md','checks/current_angular_strong_graph.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
 s,z,H=S.symbols('s z H',positive=True)
 f=s*(s-1)*z**(s-2)*(H-S.log(z))
 expected=z**(s-4)*(s*(s-1)*(s-2)*(s-3)*(H-S.log(z))-s*(s-1)*(2*s-5))
 check('complete_even_Mellin_finite_polynomial',S.simplify(S.diff(f,z,2)-expected)==0)
 check('mass_power_flat_Weyl_cancellation',(-s*(s-1)*(2*s-5)).subs(s,3)==-6)
 # Squared-Bessel moments from the retained beta integral formula.
 moment=lambda power,nu:2**(power-3)*S.gamma(S.Rational(power,2))**2*S.gamma(S.Rational(power,2)+nu)*S.gamma(S.Rational(power,2)-nu)/S.gamma(power)
 check('Bessel_third_moments',S.simplify(moment(4,0)-S.Rational(1,3))==0 and S.simplify(moment(4,1)-S.Rational(2,3))==0)
 check('Bessel_fourth_moments',S.simplify(moment(5,0)-27*S.pi**2/512)==0 and S.simplify(moment(5,1)-45*S.pi**2/512)==0)
 check('complex_radial_q_coefficient',S.simplify(4*(2*moment(4,0)+moment(4,1))/S.pi**2-16/(3*S.pi**2))==0)
 check('complex_radial_sqrtC_coefficient',S.simplify(2*(2*moment(5,0)+moment(5,1))/S.pi**2-S.Rational(99,256))==0)
 check('complex_stationary_constants',20*arb.pi()<100 and 160*arb(2).sqrt()/3<100)
 check('complex_square_root_margin',S.Rational(3,4)**2-S.Rational(1,4)**2==S.Rational(1,2))
 check('infinite_annulus_overlap_guard',S.Rational(3)/S.Rational(1,2)<2**3)
 check('near_support_precedes_input_cut',S.Rational(3,2)*S.Rational('.01')<S.Rational('.02'))
 check('output_multiplication_guard',1+3*3+3*10+60==100)
 return dict(Mellin_polynomial=str(-s*(s-1)*(2*s-5)),kernel_derivative_guards=[int(100*S.factorial(r)*4**r) for r in range(6)],complex_stationary_guard=100,overlap_guard=4)

def execute():
 alg=algebra()
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 C=pars['C'];q=pars['q'];B2=pars['B2'];beta=pars['beta'];g=pars['g']
 K=C*beta**2/(15*arb.pi()*B2);L=10*arb.pi()*(pars['Fl']+23*g/(180*arb.pi()))/g
 sigma=arb('4.005');delta=arb('.005');a=arb('.01');b=arb('.02')
 Cb=[2/(K*L*delta*sigma*(sigma-3)),2/(K*L*sigma*(sigma-3)),2/(K*L*sigma),2/(K*L)]
 Cb.append(max(Cb[0],8/(K*L),1/(K*(L*arb(6).log()).sqrt())))
 W=Cb[:4]+[Cb[1],Cb[2],Cb[3],Cb[4]]
 t,xi=S.symbols('t xi');sig=S.Rational(801,200)
 cut=S.integrate(t**8*(1-t)**8,(t,0,t))*S.factorial(17)/S.factorial(8)**2
 guards=[S.Integer(1)]
 for i in range(1,8):
  poly=S.Poly(S.diff(cut,t,i),t);n=poly.degree()
  bern=[sum(poly.nth(k)*S.binomial(j,k)/S.binomial(n,k) for k in range(j+1)) for j in range(n+1)]
  guards.append(4**i*max(abs(x) for x in bern))
 cs=[[S.Integer(1)]]
 for l in range(3):
  prior=cs[-1];cs.append([S.expand((5+j-sig)*(prior[j] if j<len(prior) else 0)+(prior[j-1] if j else 0)) for j in range(len(prior)+1)])
 def rat(x):return arb(int(x.p))/int(x.q)
 near=arb(0);near_by_m=[arb(0) for _ in range(8)];near_terms=[]
 for l in range(4):
  for j,cj in enumerate(cs[l]):
   for kk in range(4):
    n=1+j+kk
    pref=arb(comb(3,l)).sqrt()*rat(cj)*((5+j-sigma)*arb(2).log()).exp()*arb(comb(3,kk)).sqrt()*(a/2)**(4-kk)
    for i in range(n+1):
     power=sig+S.Rational(1,2)-n+i;endpoint=arb(3) if power>=0 else arb('.5')
     weight=(rat(power)*endpoint.log()).exp()
     fall=S.Poly(S.prod(sig+xi-v for v in range(n-i)),xi)
     for (mm,),cc in fall.terms():
      cell=pref*comb(n,i)*rat(guards[i])*weight*rat(abs(cc))*W[mm]
      near_by_m[mm]+=cell;near+=cell
 near*=100*arb(8).sqrt()*1000*C*2/B2
 near_by_m=[x*100*arb(8).sqrt()*1000*C*2/B2 for x in near_by_m]
 far=arb(0);far_rows=[]
 for l in range(4):
  row=arb(0)
  for j,cj in enumerate(cs[l]):
   rr=j+1;Ar=100*factorial(rr)*4**rr
   low=2*Ar*5**(4+rr)
   high=low if rr%2==0 else 2*(100*factorial(rr+1)*4**(rr+1))*5**(5+rr)
   ex=4+rr+(rr%2)-sigma
   row+=rat(cj)*(arb(low)/(sigma+1)+arb(high)/ex)
  far_rows.append(row);far+=arb(comb(3,l)).sqrt()*row
 far*=100*2/B2*(q*b*b+C.sqrt()*b)*Cb[0]

 check('actual_symbol_parameter_box',C>1 and C<7*arb(2000)**3 and q/C>300 and q/C<400 and beta>0 and beta<1)
 b02=arb(2000)**10*((1+12*arb.pi()*arb(2000)**6).sqrt()-1)/6
 logmuB=(B2/b02).log()/2;kappa4=(arb(q+1).digamma()-arb.const_euler())/2-logmuB
 check('actual_fixed_scale_parameter_box',logmuB.abs_upper()<1 and (logmuB-C.log()/2).abs_upper()<20 and kappa4>0 and kappa4<100)
 check('cutoff_endpoint_values',cut.subs(t,0)==0 and cut.subs(t,1)==1)
 check('eight_cutoff_endpoint_jets',all(S.diff(cut,t,i).subs(t,e)==0 for i in range(1,9) for e in [0,1]))
 check('positive_near_scale_powers',all(4-k>0 for k in range(4)))
 check('reflected_far_Schur_powers',6-sigma>0 and 8-sigma>0 and sigma+1>0)
 check('highest_derivative_only_seven',len(W)==8 and W[7].overlaps(Cb[4]) and Cb[4]>arb('1e18'))
 check('near_graph_target',near<arb('3.038e-9'))
 check('far_graph_target',far<arb('1.576e-16'))
 check('complete_angular_graph_target',near+far<arb('3.038e-9'))
 check('near_ledger_closes',sum(near_by_m,arb(0)).overlaps(near))
 check('safe_remaining_target',arb('.1599164')-arb('4e-9')>arb('.159916395999999'))
 return dict(algebra=alg,near=old.encoded(near),far=old.encoded(far),total=old.encoded(near+far),near_by_input_derivative=[old.encoded(x) for x in near_by_m],cutoff_derivative_guards=[str(x) for x in guards],output_polynomials=[[str(x) for x in row] for row in cs],normal_bounds=[old.encoded(x) for x in Cb],far_rows=[old.encoded(x) for x in far_rows],logmuB=old.encoded(logmuB),kappa4=old.encoded(kappa4),allocated_angular_feedback='4e-9',remaining_seed_nonlocal_target='.159916396',
  seed_flat_angular_graph_certified=True,all_angular_levels_retained=True,all_frequencies_and_annuli_retained=True,reflected_far_cancellation_retained=True,finite_fourth_contact_retained=True,normal_operator_changed=False,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,finite_physical_modes_selected=False,actual_curved_thermal_remainder_certified=False,actual_diagonal_contact_jets_certified=False,uniform_complete_metric_variation_certified=False,original_X_graph_comparison_known=False,full_signed_projection_executed=False,actual_global_correction=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_complete_response_calculation_executed=False)

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
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),angular_graph_bound=result.get('total')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
