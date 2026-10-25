"""Arithmetic supporting an ordinary frozen-response attempt; not a full certificate."""
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
STEM='current-frozen-response-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_frozen_response_feasibility.py','evidence/'+STEM+'-review.json']
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
 @lru_cache(None)
 def part(k):
  out=arb(0);by=[arb(0) for _ in range(9)]
  for n in range(k+1):
   for i in range(n+1):
    for (p,m),cc in PN[n-i].items():
     power=sig+k-n+p;geo=arb(2)**(-power);sp=H**power/(1-geo);sl=sp*((1/H).log()+arb(2).log()*geo/(1-geo))
     mm0,mm1,mm4=bounds(i,-n+i+p,m==8)
     if m==8:aa=512*mm0*WL[8]+32*(mm0+mm4)*W[7]
     else:aa=725*mm0*WL[m]+9*mm1*W[m]
     term=arb(comb(k,n)).sqrt()*comb(n,i)*abs(cc)*(aa*sp+2*mm0*W[m]*sl)
     out+=term;by[m]+=term
  return out,by
 near=arb(0);by=[arb(0) for _ in range(9)];byk={}
 for ell in range(4):
  for r in range(ell+2):
   k=4+r
   co=arb('2.02') if r==0 else arb('4.04')*cs[ell][r-1]
   outer=arb(comb(3,ell)).sqrt()*co*max(arb(1),arb(2)**(k-sig))*old.decode(series['3_'+str(k)]['upper'])
   value,bb=part(k);near+=outer*value;byk[k]=byk.get(k,arb(0))+outer*value
   for m in range(9):by[m]+=outer*bb[m]
 pref=2*C*beta**2/B2*100*4*arb(2).sqrt()*V
 ans=dict(near=old.encoded(near*pref),by_derivative=[old.encoded(v*pref) for v in by],by_k={str(k):old.encoded(v*pref) for k,v in byk.items()},V=old.encoded(V),PN_sizes=[len(x) for x in PN],scope='HIGH INPUT ALLOCATION ONLY; background-high and contacts/returns not bounded')
 def Bk(k):
  bb=arb(comb(k+3,4))/(k+2);cc=arb(comb(k+5,6))/(k+4)
  return ((arb(4*k)/(k-1)+4*(k+1)+16*bb+48*cc)/arb(2).log()+32*bb)/arb.pi()
 def DnB(n):
  out=arb(0)
  for m in range(1,3):
   for j in range(m+1):
    k=n+m-j
    if k<max(2,m):continue
    out+=([0,1,k-1][m])*comb(m,j)*Bk(k)
  return out
 R=arb('.000251');N=16
 tail=3000*arb(N+9)**14*R**N/(1-R*(arb(N+10)/(N+9))**14)
 Js={l:sum((n*DnB(n)*arb(n+1)**(l-1)*(1+arb(n+1).log()/arb(2).log())*R**(n-1) for n in range(1,N+1)),arb(0))+tail for l in range(3,8)}
 Xi=[max(rat(max(abs(a),abs(b))) for a,b in controls[i])*4**i for i in range(10)]
 En=[sum((comb(n,i)*Xi[i]*factorial(n-i)*4**(n-i) for i in range(n+1)),arb(0)) for n in range(10)]
 seed_W1=arb('7.5').sqrt()*(En[1]*En[2]).sqrt()
 @lru_cache(None)
 def geom(power):
  power=arb(power)
  geo=arb(2)**(-power);sp=H**power/(1-geo)
  return sp,sp*((1/H).log()+arb(2).log()*geo/(1-geo))
 @lru_cache(None)
 def seedpart(l):
  out=arb(0);by=[arb(0) for _ in range(8)]
  for n in range(l+1):
   for i in range(n+1):
    for (p,m),cc in PN[n-i].items():
     m0,m1,_=bounds(i,-n+i+p)
     aa=725*m0*WL[m]+9*m1*W[m]
     for amp,degree in [(140,2),(arb('1e-8'),1)]:
      sp,sl=geom(l+1-n+p+degree)
      term=amp*arb(comb(l,n)).sqrt()*comb(n,i)*abs(cc)*(aa*sp+2*m0*W[m]*sl)
      out+=term;by[m]+=term
  return 2*seed_W1*out,[2*seed_W1*v for v in by]
 seedhi=arb(0);seedbg=arb(0);seedby=[arb(0) for _ in range(8)]
 BG=[arb(0),arb(0)]
 for ell in range(4):
  for r in range(ell+2):
   k=4+r;l=k-1
   co=arb('2.02') if r==0 else arb('4.04')*cs[ell][r-1]
   outer=arb(comb(3,ell)).sqrt()*co*max(arb(1),arb(2)**(k-sig))
   val,bb=seedpart(l);seedhi+=outer*Js[l]*val
   for m in range(8):seedby[m]+=outer*Js[l]*bb[m]
   aa=sum((arb(comb(k,n)).sqrt()*2*(En[n]+En[n+1]) for n in range(k+1)),arb(0))
   bb=sum((arb(comb(k,n)).sqrt()*En[n] for n in range(k+1)),arb(0))
   for amp,powr in [(140,3),(arb('1e-8'),2)]:
    sp,sl=geom(powr)
    seedbg+=outer*Js[l]*arb(2).sqrt()*V*arb('7.5').sqrt()*amp*(aa*sp+bb*sl)
   # Bcap is sup h^(k-5/2) H_k(f_j)/(1+log(1/h)), both reflected sides.
   # The earlier unweighted sup is infinite even for the nonzero quadratic seed.
   sp,sl=geom('6.005')
   BG[0]+=outer*old.decode(series['3_'+str(k)]['derivative_upper'])*2*V**2*(sp+sl)
 p0=2*C*beta**2/B2*100
 ans2=dict(seed_gradient_high_input=old.encoded(seedhi*p0),seed_gradient_high_background=old.encoded(seedbg*p0),seed_gradient_total=old.encoded((seedhi+seedbg)*p0),seed_by_derivative=[old.encoded(v*p0) for v in seedby],gradient_series={str(k):old.encoded(v) for k,v in Js.items()},gradient_order_tail=old.encoded(tail),Bcap_coefficient=old.encoded(BG[0]*p0),Bcap_equal_half_budget_ceiling=old.encoded((arb('7.45045264695e9')-arb('521282000'))/(BG[0]*p0)),scope='cap-minus-frozen subtracted seed only; complete frozen symbol, common restorations, cap returns still present and unbounded')

 def polyadd(a,b):
  return [(a[i] if i<len(a) else arb(0))+(b[i] if i<len(b) else arb(0)) for i in range(max(len(a),len(b)))]
 def polymul(a,b):
  out=[arb(0)]*(len(a)+len(b)-1)
  for i,v in enumerate(a):
   for j,w in enumerate(b):out[i+j]+=v*w
  return out
 def scale(a,v):return [v*x for x in a]
 def DD(a):return [i*v for i,v in enumerate(a)]
 def dz(a):return [i*a[i] for i in range(1,len(a))]
 def normpoly(a,z):return sum((abs(v)*z**i for i,v in enumerate(a)),arb(0))
 delta=normpoly([arb(0)]+q[1:],arb('4.004')*H)
 check('complex_seed_inverse_contains_full_annulus',delta<arb('.001') and arb('4.004')*(1-delta/(1-delta))>4)
 check('complex_seed_mass_amplitude',abs(q[2])*arb('4.004')**2/(1+(1-delta).sqrt())<arb('1e-8') and sum((abs(q[k])*arb('4.004')**k*H**(k-3) for k in range(3,7)),arb(0))/(1+(1-delta).sqrt())<140)
 Z=q[:];jetguards=[]
 for n in range(1,5):
  jetguards.append(arb(3)**(n-1)*H**(n-1)*normpoly(Z,arb('3.003')*H)/arb('.999'))
  Z=polymul(q,dz(Z))
 check('four_normalized_actual_coordinate_jets',all(v<arb('1.01') for v in jetguards))
 check('actual_low_Fourier_guard',2*arb('1.001')**(sig+2)/arb(3).sqrt()<2)
 check('all_coordinate_power_exponents_covered',arb('1.001')**61*arb('.999')**60>1 and all(abs(sig-n+i+p)<60 for n in range(9) for i in range(n+1) for p,m in PN[n-i]))
 check('logarithmic_composition_ordinary_constants',4+384*arb(3).log()*arb(2).log()+192*arb(3).log()**2<529 and 48*(1+arb(3).log()/arb(2).log())<529 and 19*529*8<512**2 and 2*19*529*8<725**2 and 19*4<9**2)
 check('top_derivative_Hminus1_guard',(arb('3.01')+arb('1.01'))*2<16)
 check('inverse_coordinate_weak_derivative_guard',1+arb('3.003')*H*normpoly(dz(q),arb('3.003')*H)/arb('.999')<arb('1.01'))
 # Rational output coefficients share denominator q^(2+ell).
 outrows=[{0:scale(DD(q),-1),1:[arb(1)]}]
 for ell in range(3):
  row={}
  for r in range(ell+3):
   a=outrows[-1].get(r,[arb(0)])
   row[r]=polyadd(polyadd(polymul(q,polyadd(scale(a,4+r-sig),DD(a))),scale(polymul(DD(q),a),-(2+ell))),outrows[-1].get(r-1,[arb(0)]))
  outrows.append(row)
 outbounds=[]
 for ell,row in enumerate(outrows):
  vals=[]
  for r,a in row.items():
   value=normpoly(a,arb('2.002')*H)/arb('.999')**(2+ell)
   target=arb(1) if r==0 else 2*cs[ell][r-1]
   assert value<target
   vals.append(old.encoded(value))
  outbounds.append(vals)
 check('all_actual_density_and_output_returns',True)
 check('complete_gradient_series_tail',tail<arb('2.775e-35') and R*(arb(26)/25)**14<arb('.001') and Js[7]<arb('62256.431'))
 check('seed_quadratic_log_norm_obstruction',q[2]>0)
 check('executed_signed_seed_remainder',(seedhi+seedbg)*p0<arb('.000142'))
 check('executed_nonlinear_high_input_allocation',old.decode(ans['near'])<521282000)
 check('explicit_logarithmic_background_coefficient',BG[0]*p0<arb('8.123e-26'))
 budget=arb('.1599164')-arb('4.001e-9')-arb('.007452')-arb('1e-11')*881016306
 budget_B=budget-arb('8.123e-37')*arb('1e33')
 check('unchanged_physical_budget_repartition',budget.overlaps(arb('.143654232939')) and budget_B.overlaps(arb('.142841932939')))
 check('preserved_original_graph_constants',W[7]<arb('1.206e18') and WL[8]>arb('6.7402e35') and WL[8]<arb('6.7403e35'))
 check('known_allocation_not_complete_physical_certificate',True)
 return dict(outcome='unresolved_complete_strong_estimate',cap=ans,seed=ans2,
  coordinate_jet_guards=[old.encoded(v) for v in jetguards],output_coefficient_bounds=outbounds,
  actual_complete_bound_executed=False,complete_strong_estimate_status='unresolved',
  seed_cap_minus_frozen_upper='.000142',nonlinear_cap_high_input_upper=521282000,
  background_coefficient_upper='8.123e-26',background_seminorm_known=False,
  obsolete_unweighted_background_seminorm_infinite_at_actual_seed=True,
  complete_eta_star_known=False,complete_C_star_known=False,
  remaining_budget=old.encoded(budget),example_B_cap_target='1e33',remaining_with_example_B_cap=old.encoded(budget_B),
  equal_split_eta_star=old.encoded(budget_B/2),equal_split_C_star=old.encoded(budget_B/arb('2e-11')),
  best_case_B_cap_ceiling=old.encoded(budget/arb('8.123e-37')),
  optimistic_RX_alpha_squared_ceiling=old.encoded(arb('5e7')/881016306),
  optimistic_RX_alpha_ceiling=old.encoded(arb('.25')/arb('.007452')),
  all_background_orders_included=True,angular_spectrum_truncated=False,
  thermal_and_actual_returns_retained=True,new_cap_returns_certified=False,
  common_local_restorations_retained=True,complete_strong_regulator_bound_proved=False,
  original_X_graph_comparison_known=False,global_admission_known=False,
  physical_domain_changed=False,normalization_changed=False,finite_physical_modes_selected=False,
  old_unsigned_seed_estimator_failure_preserved=True,construction_obstruction=False,
  actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,
  physical_matching=False,validity_noise=False,all_sector_stability=False,
  selected_projection_screen_retained=True,projection_reevaluated=False,audit_status='closed',goal_status='active',
  next_completed_frozen_response_calculation_executed=False)

def main():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 prev=json.loads((ROOT/'evidence/current-conormal-common-regulator-attempt.json').read_bytes())['runs'][-1]
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
 print(json.dumps(dict(status=run['status'],checks=len(CHECKS),inputs=len(before),exception=run.get('exception'),seed_upper=result.get('seed_cap_minus_frozen_upper'),nonlinear_cap=result.get('nonlinear_cap_high_input_upper'),remaining=result.get('remaining_budget')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
