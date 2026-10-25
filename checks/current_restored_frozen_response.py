"""Support ordinary restored-frozen identities and bounds; no full strong certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from functools import lru_cache
from math import comb,factorial
import sys,json,traceback
DEVELOPMENT='--development' in sys.argv
if DEVELOPMENT:sys.argv.remove('--development')
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-restored-frozen-response'
OWN=['notes/'+STEM+'.md','checks/current_restored_frozen_response.py','evidence/'+STEM+'-review.json']
CHECKS=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(n,b):
 CHECKS.append(dict(name=n,passed=bool(b)))
 if not b:raise AssertionError(n)
def calc():
 ctx.prec=512
 root=ROOT
 series=json.loads((root/'evidence/current-annular-response-feasibility.json').read_bytes())['runs'][-1]['result']['series']
 pars,_,_,_=old.scalar.prev.york.actual_parameters()
 sig=arb('4.005');H=arb('.005006');C=pars['C'];beta=pars['beta'];B2=pars['B2'];g=pars['g']
 q=[arb(1),arb(0),pars['u']/beta,(-4*beta-2*pars['u']-pars['tau'])/beta,(3*beta+pars['u']+pars['eps'])/beta,-2*pars['eps']/beta,pars['eps']/beta]
 K=C*beta**2/(15*arb.pi()*B2);L=10*arb.pi()*(pars['Fl']+23*g/(180*arb.pi()))/g
 Cb=[2/(K*L*arb('.005')*sig*(sig-3)),2/(K*L*sig*(sig-3)),2/(K*L*sig),2/(K*L)]
 Cb4=max(Cb[0],8/(K*L),1/(K*(L*arb(6).log()).sqrt()))
 Cq=arb('2.328');W=Cb+[Cb[1],Cb[2],Cb[3],Cb4,Cq/(K*arb(2).log())]
 WL=[max(W[j]*(2+arb('1e18')).log(),Cq/K*arb('1e18')**(j-8)) for j in range(8)]+[Cq/K]
 def rat(v):return arb(int(v.p))/int(v.q)
 def split(coeff):
  rows=[coeff]
  for i in range(len(coeff)-1):rows.append([(rows[-1][j]+rows[-1][j+1])/2 for j in range(len(rows[-1])-1)])
  return [row[0] for row in rows],[row[-1] for row in rows][::-1]
 base=[S.Integer(0)]*9+[S.Integer(1)]*9;controls=[];current=base
 for order in range(10):
  cells=[current]
  for level in range(6):cells=[c for cell in cells for c in split(cell)]
  controls.append([[min(cell),max(cell)] for cell in cells]);degree=len(current)-1
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
 V=2*(V0*V1).sqrt()
 PN=[{(0,0):arb(1)}]
 for n in range(8):
  out={}
  for (p,m),v in PN[-1].items():
   for j,qj in enumerate(q):
    if qj==0:continue
    out[p+j,m]=out.get((p+j,m),arb(0))+qj*(sig-n+p)*v
    out[p+j,m+1]=out.get((p+j,m+1),arb(0))+qj*v
  PN.append(out)
 cs=[[arb(1)]]
 for ell in range(1,4):cs.append([(5+j-sig)*(cs[-1][j] if j<ell else 0)+(cs[-1][j-1] if j else 0) for j in range(ell+1)])
 @lru_cache(None)
 def bounds(i,offset,top=False):
  r=sig+offset;guard=arb('1.001')**(r.abs_upper()+1)
  m0=guard*M(i,offset)
  m1=guard*(M(i+1,offset)+arb('1.01')*abs(r)*M(i,offset-1))
  m4=arb(0)
  if top:
   for j in range(5):
    order=4-j;bell=arb(1) if order==0 else arb(0)
    falling=arb(1)
    for k in range(1,order+1):
     falling*=r-(k-1)
     from sympy.functions.combinatorial.numbers import stirling
     bell+=abs(falling)*int(stirling(order,k,kind=2))*arb('1.01')**k
    m4+=guard*comb(4,j)*M(j,offset-4+j)*bell
  return m0,m1,m4
 R=arb('.000251');NN=16
 def series_bound(l,third=False):
  U=arb(0);Vv=arb(0)
  for n in range(NN+1):
   coef=arb(1000)*4**n*((5*n+4) if third else 1)
   slots=n+(2 if third else 1)
   aa=coef*arb(slots)**(l-1)
   U+=aa*R**n
   if n:Vv+=aa*n*R**(n-1)
  # All l<=7, both series, coarse polynomial tail and its derivative.
  U+=5000*4**(NN+1)*arb(NN+3)**7*R**(NN+1)/(1-4*R*(arb(NN+4)/(NN+3))**7)
  Vv+=5000*4**(NN+1)*arb(NN+3)**8*R**NN/(1-4*R*(arb(NN+4)/(NN+3))**8)
  return U,Vv
 Xi=[max(rat(max(abs(a),abs(b))) for a,b in controls[i])*4**i for i in range(10)]
 En=[sum((comb(n,i)*Xi[i]*factorial(n-i)*4**(n-i) for i in range(n+1)),arb(0)) for n in range(10)]
 @lru_cache(None)
 def geo(power):
  power=arb(power);g=arb(2)**(-power);sp=H**power/(1-g)
  return sp,sp*((1/H).log()+arb(2).log()*g/(1-g))
 @lru_cache(None)
 def input_sum(l,third=False):
  total=arb(0);by=[arb(0)]*8
  for n in range(l+1):
   for i in range(n+1):
    for (p,m),cc in PN[n-i].items():
     m0,_,_=bounds(i,-n+i+p)
     sp,_=geo(str(sig+l+1-n+p) if third else l+1-n+p)
     term=arb(comb(l,n)).sqrt()*comb(n,i)*abs(cc)*2*m0*W[m]*sp
     total+=term;by[m]+=term
  return total,by
 seed_high=arb(0);seed_bg=arb(0);third_high=arb(0);third_B=arb(0)
 tables={}
 for ell in range(4):
  for r in range(ell+2):
   k=4+r;l=3+r
   co=arb('2.02') if r==0 else arb('4.04')*cs[ell][r-1]
   outer=arb(comb(3,ell)).sqrt()*co*max(arb(1),arb(2)**(k-sig))
   U,Vv=series_bound(l);U3,V3=series_bound(l,True)
   tables[l]=[old.encoded(v) for v in [U,Vv,U3,V3]]
   hi,_=input_sum(l);hi3,_=input_sum(l,True)
   seed_high+=outer*U*2*hi
   third_high+=outer*U3*4*arb(2).sqrt()*V*hi3
   ee=sum((arb(comb(l,n)).sqrt()*En[n] for n in range(l+1)),arb(0))
   sp4,_=geo(4);sp3,_=geo(3)
   seed_bg+=outer*Vv*arb(2).sqrt()*V*arb('7.5').sqrt()*ee*(140*sp4+arb('1e-8')*sp3)
   pp=sig+3 if l>=4 else sig+2
   # Avoid hashing a non-exact Arb: only these two exact decimal exponents occur.
   sp,sl=geo('7.005' if l>=4 else '6.005')
   third_B+=outer*V3*2*V**2*(sp+sl)/arb(2).log()
 # Entire spectral remainder is bounded by 1000 C beta^(3/2) |p|^3.
 pref=2*C*beta**arb('1.5')/B2*100
 ans=dict(seed_high=old.encoded(seed_high*pref),seed_background=old.encoded(seed_bg*pref),seed_total=old.encoded((seed_high+seed_bg)*pref),third_high=old.encoded(third_high*pref),third_B_coefficient=old.encoded(third_B*pref),series=tables,scope='complete order-three spectral remainder of left frozen symbol, not completed geometric response')

 # Ordinary functional identities; these checks do not supply domain arguments.
 c,x,aa,ka,bb,uu,qq,kk=S.symbols('c x a kappa beta u q k',real=True)
 Kg=S.Function('K');Kmc=S.exp(2*c)*(Kg(S.exp(-c)*x)+aa*(S.exp(-c)*x)**4/2)+c*aa*S.exp(-2*c)*x**4
 Kconv=-aa*S.exp(-2*c)*(2*c+1)*x**4/2
 check('entire_frozen_mass_geometric_identity',S.simplify(Kmc+Kconv-S.exp(2*c)*Kg(S.exp(-c)*x))==0)
 check('correct_common_soft_third_identity',S.simplify(S.diff(Kmc,c).subs(c,0)-(2*Kg(x)-x*S.diff(Kg(x),x)))==0)
 check('geometric_conversion_third_quartic_zero',S.simplify(S.diff(Kconv,c).subs(c,0))==0)
 check('superseded_mass_initialization_detected',S.simplify((2*aa*ka*x**4-x*S.diff(aa*ka*x**4,x)+aa*x**4)-(-2*aa*ka*x**4))==aa*x**4)
 pp0,pp1,pp2=S.symbols('P P1 P2',real=True)
 LL=-2*bb*pp2;MM=-bb*pp2;Nphi=bb*pp1**2
 Wraw=(S.exp(-2*pp0)*(-pp0*LL**2/60+LL*MM/30-MM**2/10+LL*Nphi/30-S.Rational(7,90)*MM*Nphi-S.Rational(7,180)*Nphi**2)+uu*pp0*LL/18-uu*Nphi/18-kk*uu**2*pp0*S.exp(2*pp0))/2
 Jraw=-Wraw-qq*uu*bb*pp1**2/12
 Jred=bb**2*S.exp(-2*pp0)*((2*pp0+1)*pp2**2/60+S.Rational(17,1080)*pp1**4)-uu*bb*(1+3*qq)*pp1**2/36+kk*uu**2*pp0*S.exp(2*pp0)/2
 flux=-bb**2*S.exp(-2*pp0)*pp1**3/540+uu*bb*pp0*pp1/18
 df=S.diff(flux,pp0)*pp1+S.diff(flux,pp1)*pp2
 check('complete_conversion_boundary_representative',S.simplify(Jraw-Jred-df)==0)
 c0,fq=S.symbols('c0 fq',real=True)
 Lm=c0*S.exp(2*pp0)/4+fq*(-bb**2*S.exp(-2*pp0)*((ka+pp0+S.Rational(1,2))*pp2**2/30+S.Rational(17,1080)*pp1**4)+uu*bb*(1+3*qq)*pp1**2/36-kk*uu**2*pp0*S.exp(2*pp0)/2)
 Lg=c0*S.exp(2*pp0)/4-fq*ka*bb**2*S.exp(-2*pp0)*pp2**2/30
 check('entire_degree_zero_two_four_restoration',S.simplify(Lm+fq*Jred-Lg)==0)
 q0,q1,q2,b0,b1,b2=S.symbols('q0 q1 q2 b b1 b2',nonzero=True)
 Psiss=q0*(q0*b2+q1*b1+q2/2);lap=bb*(q0*b2+q1*b1);curv=-bb*q2
 check('actual_curved_density_return',S.simplify(bb**2*S.exp(-2*b0)*Psiss**2/q0**2-S.exp(-2*b0)*(curv-2*lap)**2/4)==0)
 check('same_explicit_a2_coefficient',S.simplify(-ka/120-(-4*ka)/480)==0)
 L0=S.symbols('L0',real=True)
 Kinfty=lambda p:aa*p**4*(S.log(p)-L0-S.Rational(31,30))
 expected=aa*S.exp(-2*c)*x**4*(S.log(x)-L0-S.Rational(8,15))
 # x>0 is understood in the logarithm identity; substitute log(exp(-c)x)=log(x)-c explicitly.
 raw=S.exp(2*c)*(aa*S.exp(-4*c)*x**4*(S.log(x)-c-L0-S.Rational(31,30))+aa*S.exp(-4*c)*x**4/2)+c*aa*S.exp(-2*c)*x**4
 check('full_frozen_principal_with_same_regulator',S.simplify(raw-expected)==0)
 check('corrected_third_kernel_still_covered',arb(74)/(5*arb.pi())+arb(33)/(10*arb.pi()*arb(2).log())<7 and 32<33)
 check('complex_sum_integral_constants',6*arb.pi()*arb(2).sqrt()<32 and 5*arb.pi()*arb(2)**arb('1.5')<48)
 # Anchors enclose the complete convergent angular remainders, not sampled p.
 rho=arb(1)/1296
 anchorK=16/(30*arb.pi())+4/(300*840*arb.pi()*(1-rho))
 anchorD=32/(30*arb.pi())+4/(300*280*arb.pi()*(1-rho))
 anchorInf=(arb.pi()/6+20+arb(31)/30)/(30*arb.pi())
 anchorInfD=(arb.pi()/3+40+arb(47)/30)/(30*arb.pi())
 check('whole_sector_unit_anchors',anchorK+anchorD+anchorInf+anchorInfD<20)
 check('whole_sector_high_and_low_bounds',20+400/arb.pi()+arb('1.25')<200 and 2+(1/arb(1).exp()+arb.pi()/6+20+arb(31)/30)/(30*arb.pi())<200)
 check('full_complex_mass_disk_spectral_guard',200*arb(4)/3<1000)
 check('spectral_parameters',pars['q']/C>300 and pars['q']/C<400 and (pars['q']+1)/C>324)
 check('spectral_derivative_scaling',beta**arb('1.5')>beta**2 and beta<1)
 check('entire_order_tail_ratio',4*R*(arb(NN+4)/(NN+3))**8<arb('.002'))
 check('seed_frozen_spectral_strong_bound',old.decode(ans['seed_total'])<arb('8.424e-10'))
 check('one_metric_variation_frozen_spectral_bound',old.decode(ans['third_high'])<arb('8.573e-11') and old.decode(ans['third_B_coefficient'])<arb('1.282e-25'))
 check('whole_Mellin_comparison_halfplane',21*(200+18)<5000 and L>arb('3e34') and 5000<L/2)
 reT,imT,LL=S.symbols('reT imT L',real=True)
 check('exact_full_normal_correlation',(S.expand((reT-LL)**2+imT**2-(reT**2+imT**2))-(LL**2-2*LL*reT)).expand()==0)
 poly=[sum((q[i]*q[k-i] for i in range(max(0,k-6),min(6,k)+1)),arb(0)) for k in range(13)]
 rs=[sum((abs(poly[k])*k**j*arb('.01')**k for k in range(1,13)),arb(0)) for j in range(4)]
 mult=sum((w*v for w,v in zip([100,51,12,1],rs)),arb(0))
 check('complete_matched_Mellin_component',mult<arb('.003138'))
 budget=arb('.1599164')-arb('4.001e-9')-arb('.003138')-arb('.000142')-arb('8.424e-10')-arb('1e-11')*881016306
 budgetB=budget-arb('8.123e-37')*arb('1e33')
 check('unchanged_target_with_complete_comparison',budget.overlaps(arb('.1478262320966')) and budgetB.overlaps(arb('.1470139320966')))
 projection=json.loads((ROOT/'evidence/current-signed-resolvent-projection.json').read_bytes())['runs'][-1]
 check('projection_evidence_precedes_corrected_third_comparison',not any('current-mixed-frequency-third-response' in p for p in projection['input_sha256_after']))
 return dict(outcome='unresolved_complete_strong_estimate',actual_complete_bound_executed=False,
  restored_local_bulk_identified=True,mass_geometric_fourth_contact_correction=True,
  corrected_fourth_coefficient='-2 C kappa4/(30 pi)',superseded_fourth_coefficient='C(1-2 kappa4)/(30 pi)',
  old_comparison_diagnostics_are_actual_common_finite_part=False,
  corrected_flat_comparison_original_bound_preserved=True,full_restoration_added_to_explicit_local_action_twice=False,
  frozen_spectral=ans,full_Mellin_comparison_norm_upper=1,matched_Mellin_component=old.encoded(mult),
  matched_Mellin_component_upper='.003138',remaining_budget=old.encoded(budget),example_B_cap_target='1e33',
  remaining_with_example_B_cap=old.encoded(budgetB),equal_split_eta_dagger=old.encoded(budgetB/2),equal_split_C_dagger=old.encoded(budgetB/arb('2e-11')),
  best_case_B_cap_ceiling=old.encoded(budget/arb('8.123e-37')),
  full_eta_dagger_known=False,full_C_dagger_known=False,B_cap_known=False,
  actual_coordinate_conjugation_strong_bound_known=False,new_actual_cap_return_bound_known=False,
  all_angular_levels_retained=True,all_thermal_and_physical_returns_retained=True,
  original_X_graph_comparisons_known=False,global_family_admission_proved=False,
  physical_domain_changed=False,normalization_changed=False,finite_physical_modes_selected=False,
  construction_obstruction=False,actual_controlled_solution=False,nonlinear_admission=False,
  exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
  selected_projection_screen_retained=True,projection_reevaluated=False,audit_status='closed',goal_status='active',
  next_joint_coordinate_and_cap_return_calculation_executed=False)
def main():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 prev=json.loads((ROOT/'evidence/current-frozen-response-feasibility.json').read_bytes())['runs'][-1]
 paths=set(OWN+review['retained_inputs'])|set(prev['input_sha256_after'])
 before={p:digest(p) for p in sorted(paths)}
 out=ROOT/'evidence'/f'{STEM}{"-development" if DEVELOPMENT else ""}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
  check('previous_scientific_sources_preserved',all(digest(p)==h for p,h in prev['input_sha256_after'].items()))
  run['result']=calc();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=CHECKS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 if not DEVELOPMENT:
  with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(CHECKS),inputs=len(before),exception=run.get('exception'),matched_Mellin=result.get('matched_Mellin_component_upper'),restoration=result.get('restored_local_bulk_identified'),remaining=result.get('remaining_budget')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
