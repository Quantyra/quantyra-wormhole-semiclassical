"""Coupled classical local graph bound, on the actual finite-seed collar."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as prev
STEM='current-coupled-local-graph-bound'
OWN=['notes/'+STEM+'.md','checks/current_coupled_local_graph_bound.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={}
def check(name,ok,detail=''):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def execute():
 z,s,t=S.symbols('z s t');be,uu,ep,ta,ell,mag,cr=S.symbols('beta u eps tau ell mag cr')
 Q=be+uu*z**2+(-4*be-2*uu-ta)*z**3+(3*be+uu+ep)*z**4-2*ep*z**5+ep*z**6
 # Direction z^s removed in these symbolic row pencils; derivation rules
 # differentiate the full product before removing it.
 der=lambda v:S.diff(v,z)+s*v/z
 lap=lambda v:der(Q*der(v))
 D0=-S.diff(Q,z,2)+6*S.diff(Q,z)/z-12*Q/z**2
 def rows(n,b):
  d=-2*lap(b)-6*lap(n)+12*Q*der(n)/z
  v=d-2*D0*b
  FnR=-12*der(der(Q*v))-der((24*Q/z-12*S.diff(Q,z))*v-12*Q*(4*uu+2*D0)*der(n))
  FbR=(16*uu**2+4*D0**2)*b-4*D0*d+4*der(S.diff(Q,z)*v)-4*der(der(Q*v))
  FnE=z**-2*(2*D0*n+d)+4*uu*z**-2*(n+b)+ell*z**-4*(48*n+24*b)+der(der(-6*Q*z**-2*n))-der((12*Q/z-6*S.diff(Q,z))*z**-2*n-6*Q*z**-2*der(n))
  FbE=4*uu*z**-2*(n+b)+ell*z**-4*(24*n+12*b)-2*der(der(Q*z**-2*n))+2*der(S.diff(Q,z)*z**-2*n)-4*mag*b
  return S.expand(z**4*(FnE+cr*FnR)),S.expand(z**4*(FbE+cr*FbR))
 rn=rows(1,0);rb=rows(0,1)
 TT=lambda row:S.expand(z*S.diff(row[1],z)+(s-4)*row[1]-(z*S.diff(row[0],z)+(s-6)*row[0])/3)
 an=rn[0];ab=rb[0];tn=TT(rn);tb=TT(rb)
 a=(s+3)/(3*(s+1))
 matrix=[[an,S.cancel(ab-a*an)],[tn,S.cancel(tb-a*tn)]]
 pp,h,zh,env=prev.scalar.prev.york.actual_parameters();B2=pp['B2'];C=pp['C'];pi=arb.pi();r=(96*pi*arb('1e22')-1/(60*pi))*pp['g'];Cv=1+11*pp['g']/(180*pi);sig=arb('4.005');delta=arb('.005');a0=Cv-r*sig*(sig-3)
 K=C*pp['beta']**2/(15*pi*B2);L=10*pi*(pp['Fl']+23*pp['g']/(180*pi))/pp['g']
 ct=[1/(6*pp['beta']*a0*delta*(sig+1)),1/(6*pp['beta']*a0*(sig+1)),1/(6*pp['beta']*a0),1/(12*pp['beta']*(a0*r).sqrt()),1/(6*pp['beta']*r)]
 cb=[2/(K*L*delta*sig*(sig-3)),2/(K*L*sig*(sig-3)),2/(K*L*sig),2/(K*L)]
 mapping={be:pp['beta'],uu:pp['u'],ep:pp['eps'],ta:pp['tau'],ell:arb(2000)**-10*B2/C,mag:pi*pp['q']**2*arb(2000)**8/(C*B2),cr:8*pi*arb('1e22')*C/B2}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(v) for v in e.args),arb(0))
  if e.is_Mul:
   v=arb(1)
   for w in e.args:v*=val(w)
   return v
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 def cutoff_constants(k):
  v=S.symbols('v');cut=1-35*v**4+84*v**5-70*v**6+20*v**7
  caps=[];l2=[]
  for j in range(4):
   poly=S.Poly(cut,v);power=[poly.nth(i) for i in range(8)];cells=[]
   for cell in range(128):
    lo=S.Rational(cell,128);width=S.Rational(1,128)
    local=[sum(power[h]*comb(h,i)*lo**(h-i)*width**i for h in range(i,8)) for i in range(8)]
    bern=[sum(local[h]*S.Rational(comb(i,h),comb(7,h)) for h in range(i+1)) for i in range(8)]
    cells.append((val(max(abs(w) for w in bern))*((arb(cell+1)/128-1)*k).exp()).abs_upper())
   trans=max(cells);left=arb(k)**j*(-arb(k)).exp();caps.append(max(trans,left.abs_upper()))
   l2.append((trans**2+left**2/(2*k)).sqrt())
   cut=S.expand(S.diff(cut,v)+k*cut)
  return sum((comb(3,j)*caps[j] for j in range(4)),arb(0)),sum((comb(3,j)*l2[j] for j in range(4)),arb(0)),caps
 cuts={k:cutoff_constants(k) for k in range(2,13)}
 check('actual_positive_scalar_denominator',a0>arb('.999999'))
 check('whole_frequency_shape_reserve',L>10000)
 check('normal_scalar_conversion', (pp['g']*B2-C*pp['beta']).contains(0))
 out=[]
 for i in range(2):
  row=[]
  for j in range(2):
   num=S.expand(matrix[i][j]*(s+1)) if j else S.expand(matrix[i][j]);num=S.cancel(num).expand();pol=S.Poly(num,z)
   check('positive_powers_only_%d_%d'%(i,j),all(k[0]>=0 for k in pol.monoms()))
   cols=[];total=arb(0);low_total=arb(0);high_total=arb(0)
   for (k,),cp in pol.terms():
    if k==0:continue
    cp=S.Poly(S.expand(cp),s)
    if j:
     quot,rem=S.div(cp,S.Poly(s+1,s));p=quot.as_expr();extra=val(rem.as_expr()).abs_upper()*cb[0]/(sig+1)
    else:p=cp.as_expr();extra=arb(0)
    # Expansion in i xi, with the actual sigma and all cancellations.
    p=S.Poly(S.expand(p.subs(s,t+S.Rational(801,200))),t)
    caps=cb if j else ct
    check('entry_order_%d_%d_%d'%(i,j,k),p.degree()<len(caps))
    bound=extra+sum((val(p.nth(v)).abs_upper()*caps[v] for v in range(p.degree()+1)),arb(0))
    mult=arb('.01')**k*(k**3+12*k**2+51*k+100)
    total+=bound*mult
    # Certified low/high split: no sampled Fourier frequencies.
    low_symbol=sum((val(p.nth(v)).abs_upper() for v in range(p.degree()+1)),arb(0))
    if j:low_symbol+=val(rem.as_expr()).abs_upper()/(sig+1)
    low_fourier=low_symbol*(cb[1] if j else ct[1])/(2*delta).sqrt()
    hi_caps=[cb[3]]*4 if j else [ct[2],ct[2],ct[2],ct[3],ct[4]]
    high_fourier=sum((val(p.nth(v)).abs_upper()*hi_caps[v] for v in range(p.degree()+1)),arb(0))
    if j:high_fourier+=val(rem.as_expr()).abs_upper()*cb[3]/(sig+1)
    mm,ll,_=cuts[k];low_total+=arb('.01')**k*ll*low_fourier;high_total+=arb('.01')**k*mm*high_fourier
    cols.append({'z_power':k,'symbol':str(cp.as_expr()),'fourier_bound':prev.encoded(bound),'localized_bound':prev.encoded(bound*mult),'low':prev.encoded(arb('.01')**k*ll*low_fourier),'high':prev.encoded(arb('.01')**k*mm*high_fourier)})
   sharp=(low_total**2+high_total**2).sqrt()
   row.append({'total':prev.encoded(total),'sharp':prev.encoded(sharp),'low':prev.encoded(low_total),'high':prev.encoded(high_total),'terms':cols})
  out.append(row)
 data={'matrix_bounds':out,'scalar_multiplier_constants':[prev.encoded(v) for v in ct],'shape_multiplier_constants':[prev.encoded(v) for v in cb],'cutoff_constants':{k:[prev.encoded(v[0]),prev.encoded(v[1])] for k,v in cuts.items()}}

 sharp=[[prev.decode(c['sharp']) for c in row] for row in out]
 coarse=[[prev.decode(c['total']) for c in row] for row in out]
 eta=max(sum(row,arb(0)).abs_upper() for row in sharp)
 check('exact_cancellation_local_block_target',eta<arb('.06'))
 check('sharper_local_guard',eta<arb('.052431'))
 check('coarse_bound_does_not_certify',sum(coarse[1],arb(0))>1)
 pr=prev.decode(json.loads((ROOT/'evidence/current-outer-graph-feedback.json').read_bytes())['runs'][-1]['result']['actual_log_principal_feedback_upper'])
 eta_remaining=arb('.187');joint_known=arb('1e-14');affine=arb('7.49e-12');total=eta+pr+eta_remaining
 conditional=(affine+total*joint_known)/(1-total)
 check('remaining_allocation_fits',eta+eta_remaining<arb('.24'))
 check('conditional_joint_return_target',conditional<arb('1e-11'))
 data.update(local_block_row_sum_upper=prev.encoded(eta),prior_quantum_principal_upper=prev.encoded(pr),
  remaining_feedback_threshold=prev.encoded(arb('.24')-eta),selected_remaining_feedback='.187',
  conditional_total_feedback=prev.encoded(total),conditional_joint_E_upper=prev.encoded(conditional),
  conditional_joint_known_source_cap='1e-14',conditional_joint_affine_cap='7.49e-12',
  whole_classical_local_block_certified=True,full_graph_feedback_certified=False,
  joint_affine_and_known_source_estimates_certified=False,shape_only_k3_substituted_for_joint_known_source=False,
  actual_E3_upper_known=False,actual_Efar_upper_known=False,actual_global_correction_computed=False,
  actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,
  semiclassical_validity_noise=False,all_sector_stability=False)
 return data

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
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['symbolic']=SYMBOLIC;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
