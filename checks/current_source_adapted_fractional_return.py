"""Actual reconstruction map endpoint constants and full fractional residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import argparse,json,sys,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_source_adapted_local_return as lr
arb,ctx,dec,enc,J=lr.arb,lr.ctx,lr.dec,lr.enc,lr.J
ROOT=lr.ROOT;STEM='current-source-adapted-fractional-return'
INPUTS=['notes/'+STEM+'.md','checks/current_source_adapted_fractional_return.py',
 'checks/current_source_adapted_local_return.py','checks/current_source_adapted_geodesic_first_attempt.py',
 'checks/current_localized_first_residual_range_test.py','checks/current_localized_short_time_trial.py',
 'checks/current_quartic_ritz_rank_test.py',
 'notes/current-source-adapted-geodesic-first-attempt.md','notes/current-source-adapted-local-return.md',
 'notes/current-geodesic-two-sided-remainder.md','notes/current-seed-full-hessian-realization.md',
 'evidence/current-source-adapted-local-return-attempt4.json',
 'evidence/current-source-adapted-local-return-verification.json',
 'evidence/current-source-adapted-geodesic-first-attempt.json',
 'evidence/current-source-adapted-geodesic-first-attempt-verification.json',
 'notes/current-first-quantum-solve-attempt.md','evidence/current-first-quantum-solve-attempt.json',
 'evidence/current-global-green-trial-enclosure.json','evidence/current-global-mass-conversion.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p):return sha256(p.read_bytes()).hexdigest()
def up(x):return arb(abs(x).upper())
def maxup(a,b):return arb(max(up(a).upper(),up(b).upper()))
def sqrtjet(x):
 out=[x.c[0].sqrt()]
 assert out[0].is_finite() and not out[0].contains(0)
 for n in range(1,len(x.c)):
  out.append((x.c[n]-sum((out[j]*out[n-j] for j in range(1,n)),arb(0)))/(2*out[0]))
 return J(out)

def load():
 p={k:dec(v) for k,v in lr.result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 p['At']=dec(lr.result('current-global-green-trial-enclosure')['constants']['A_t'])
 c=lr.constants('current-source-adapted-geodesic-first-attempt')
 mp={k:dec(v) for k,v in lr.result('current-global-mass-conversion')['parameters'].items()}
 L,a,ga,beta=[p[k] for k in ['L','a','gamma','beta']]
 u,eps,tau=[mp[k] for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+eps,-2*eps,eps]
 bs=[sum((Q[j]*comb(j,i)*L**(j-i)*(-1)**i for j in range(i,7)),arb(0)) for i in range(7)]
 assert bs[0].contains(0) and bs[1].overlaps(2*ga) and bs[2].overlaps(a)
 bs[:3]=[arb(0),2*ga,a]
 assert a>arb('.999') and a<arb('1.001') and bs[3]>arb('-1.34') and bs[3]<arb('-1.32')
 assert bs[4]>arb('.49') and bs[4]<arb('.51') and abs(bs[5])+abs(bs[6])<arb('1e-25')
 assert ga<arb('1e-30') and L<=1 and L>arb('.99') and c['muY']>1
 return p,c,bs

def geometry_jets(p,bs,x):
 L,a,ga=[p[k] for k in ['L','a','gamma']]
 z=L-ga*x;q2=2*x+a*x*x
 R=1+sum((bs[j]*ga**(j-2)*x**(j-1)/(2+a*x) for j in range(3,7)),0)
 return z,q2,R,q2*R

def scalar_bounds(F,m,n=5):
 bell=[[arb(0)]*(n+1) for _ in range(n+1)];bell[0][0]=arb(1)
 for j in range(1,n+1):
  for k in range(1,j+1):
   bell[j][k]=sum((comb(j-1,i-1)*m[i]*bell[j-i][k-1] for i in range(1,j-k+2)),arb(0))
 out=[up(F.c[0])]
 for j in range(1,n+1):
  out.append(up(sum((factorial(k)*up(F.c[k])*bell[j][k] for k in range(1,j+1)),arb(0))))
 return out

def cell_scalar(p,c,bs,s0,s1):
 a,ga,L=p['a'],p['gamma'],p['L']
 x0=s0/ga;x1=s1/ga
 x=J([arb((x0+x1)/2,(x1-x0)/2),arb(1)]+[arb(0)]*9)
 z,q2,R,q=geometry_jets(p,bs,x)
 assert R.c[0]>0
 M=z*q.d().d()+4*ga*q.d()
 A=z**5*R/L**5;B=z**4*(2*M/3-z*R*c['muY'])/L**5
 C=z**4*(5*R-z*R.d()/ga)/L**5;E=-2*z**4/(3*L**5)
 D=z*q.d().d().d()/ga+3*q.d().d()
 root=up(q2.c[0]).sqrt();b=up(1+a*x.c[0]);rt2=arb(2).sqrt()
 mx=[arb(0),root,rt2*b,rt2*a*root,2*a*b,2*a*a*root,2*rt2*a*a*b]
 return {name:scalar_bounds(F,mx,4 if name=='D' else 5) for name,F in [('A',A),('B',B),('C',C),('E',E),('D',D)]}

def cell_connection(p,bs,s0,s1):
 a,ga=p['a'],p['gamma']
 h0=(s0/ga).sqrt();h1=(s1/ga).sqrt()
 h=J([arb((h0+h1)/2,(h1-h0)/2),arb(1)]+[arb(0)]*9)
 S=sqrtjet(2+a*h*h)
 U=sum((bs[j]*ga**(j-2)*h**(2*j-6) for j in range(3,7)),0)/(2+a*h*h)
 R=1+h**4*U
 assert R.c[0]>0
 b=S*(4*U+h*U.d())/(4*R)
 cc=-U*(2+h**4*U)*(1+a*h*h)/S-R*S*(4*U+h*U.d())/4
 terms={(0,0,0):-b,(1,0,1):b,(1,1,0):b,(0,1,1):cc}
 rot={(0,1):-1,(1,0):1};out=[]
 for j in range(4):
  hp=up(h.c[0]**(3-j))
  out.append(up(sum(((hp*up(t.c[0]))**2 for t in terms.values()),arb(0)).sqrt()))
  if j==3:break
  new={}
  def add(key,val):new[key]=new[key]+val if key in new else val
  for key,t in terms.items():
   add(key+(0,),S*((3-j)*t+h*t.d())/2)
   for k in [0,1]:
    if (k,key[0]) in rot:add((k,)+key[1:]+(1,),rot[k,key[0]]*(1+a*h*h)/S*t)
   for index in range(1,len(key)):
    for k in [0,1]:
     if (key[index],k) in rot:
      kk=list(key);kk[index]=k
      add(tuple(kk)+(1,),-rot[key[index],k]*(1+a*h*h)/S*t)
  terms=new
 return out

def full_map_constants(p,sc,dc):
 a,ga,L=p['a'],p['gamma'],p['L'];rt2=arb(2).sqrt()
 root=(2*ga*L+a*L*L).sqrt();b=ga+a*L
 ms=[arb(0),root,rt2*b,rt2*a*root,2*a*b,2*a*a*root,2*rt2*a*a*b]
 dj=[arb(1),arb(1)]
 for j in range(2,7):dj.append(dj[-1]*maxup(arb(1),4*(j-1)**2*a).sqrt())
 # The resolvent multiplier is <=1 for all needed derivative orders.
 W=[L*sc['D'][0]]
 for j in range(1,6):
  val=arb(0)
  for i in range(j):
   for k in range(j-i):
    ell=j-1-i-k
    val+=comb(j-1,i)*comb(j-1-i,k)*sc['D'][i]*dj[k]*ms[ell+1]
  W.append(up(val))
 source=[]
 for j in range(6):
  v=sum((comb(j,i)*(sc['A'][i]*dj[j-i]+sc['B'][i]*dj[j-i]+sc['E'][i]*W[j-i]) for i in range(j+1)),arb(0))
  for i in range(j+1):
   for k in range(j-i+1):
    ell=j-i-k
    v+=comb(j,i)*comb(j-i,k)*sc['C'][i]*ms[k+1]*dj[ell+1]
  source.append(up(v))
 # Mixed reference/actual covariant derivatives; all connection terms retained.
 T={(0,j):source[j] for j in range(6)}
 for k in range(5):
  for j in range(5-k):
   value=T[k,j+1]
   if k:
    value+=k*sum((comb(j,i)*dc[i]*T[k,j-i] for i in range(j+1)),arb(0))
   T[k+1,j]=up(value)
 C4=100*T[4,0];C5=arb(10)**(arb(5)/2)*T[5,0]
 return C4,C5,dict(model_derivative_constants=dj,s_derivative_bounds=ms,
  W_derivative_constants=W,reconstructed_reference_derivative_constants=source,
  mixed_derivative_constants={str(k):v for k,v in T.items()})

def execute():
 ctx.prec=512;p,c,bs=load();L,ga=p['L'],p['gamma']
 edges=[arb(0),ga]
 while edges[-1]<arb('.001'):
  edges.append(arb((2*edges[-1]).mid()))
 while edges[-1]+arb('.002')<L:
  edges.append(edges[-1]+arb('.002'))
 edges.append(L)
 sc={n:[arb(0)]*(5 if n=='D' else 6) for n in ['A','B','C','E','D']}
 dc=[arb(0)]*4;cells=[];queue=list(zip(edges,edges[1:]));processed=0
 while queue:
  s0,s1=queue.pop(0)
  try:
   ss=cell_scalar(p,c,bs,s0,s1);dd=cell_connection(p,bs,s0,s1)
   assert all(v.is_finite() for row in ss.values() for v in row) and all(v.is_finite() for v in dd)
  except (AssertionError,ValueError):
   assert s1-s0>ga/1024,'Cannot certify this whole-cell geometry'
   mid=(s0+s1)/2;queue[0:0]=[(s0,mid),(mid,s1)];continue
  for n in sc:
   sc[n]=[maxup(v,w) for v,w in zip(sc[n],ss[n])]
  dc=[maxup(v,w) for v,w in zip(dc,dd)]
  cells.append(dict(s0=enc(s0),s1=enc(s1),scalar_bounds={n:[enc(v) for v in row] for n,row in ss.items()},
   connection_bounds=[enc(v) for v in dd]))
  processed+=1
  if processed%100==0:print(json.dumps(dict(cells=processed,queued=len(queue),D=[v.str(8) for v in dc])),flush=True)
 C4,C5,detail=full_map_constants(p,sc,dc);flow=(C4*C5).sqrt()
 actual=flow*c['original_model_fractional_scale'];Egeo=actual*c['actual_reference_fractional_output_coefficient']
 local=json.loads((ROOT/'evidence/current-source-adapted-local-return-attempt4.json').read_bytes())
 lc={k:dec(v) for k,v in local['result']['constants'].items()}
 source=up(c['complete_actual_source_debit']);model=up(c['complete_model_density_remainder']);ret=up(lc['complete_local_return'])
 total=source+model+ret+Egeo
 account=arb('4.308605e32');target=arb('9.955e32')
 first=lr.constants('current-first-quantum-solve-attempt')
 e0=up(first['complete_geodesic_form_error'])
 reference_energy_error=total/(1-e0)
 stage=first['conditional_common_stage_output_error_ceiling']
 trial_energy_upper=(up(first['proved_source_dual_upper'])+total)/(1-e0)
 coarse_actual_difference=e0*trial_energy_upper
 constants=dict(C4=C4,C5=C5,actual_reconstruction_fractional_comparison=flow,
  actual_fractional_flow_norm_upper=actual,complete_reference_nonlocal_response=Egeo,
  quantum_debit_over_remaining_account=Egeo/account,
  complete_exact_analytic_first_reference_residual=total,
  first_residual_over_sufficient_screen=total/target,
  remaining_for_UNEXECUTED_residual_application_errors=account-Egeo,
  reference_inverse_energy_error=reference_energy_error,
  analytic_energy_error_over_conditional_common_stage_allowance=reference_energy_error/stage,
  remaining_for_UNEXECUTED_caller_Et_errors=stage-reference_energy_error,
  reference_energy_error_bolt_coefficient_product_NOT_complete_output=c['first_residual_to_bolt_multiplier']*reference_energy_error,
  repeated_source_error_component_with_known_stage_amplification_NOT_complete_output=
   c['first_residual_to_bolt_multiplier']*first['stage_error_amplification_sum']*reference_energy_error,
  actual_trial_Et_upper=trial_energy_upper,
  known_coarse_actual_minus_reference_upper=coarse_actual_difference,
  coarse_actual_difference_upper_over_complete_first_screen=coarse_actual_difference/arb('1e33'))
 for v in constants.values():assert v.is_finite()
 status='fractional_map_enclosed_exact_trial_screen_passes_numerical_applications_open' if total<target and Egeo<account else 'fractional_map_estimator_excessive'
 def encode_detail(x):
  if isinstance(x,dict):return {k:encode_detail(v) for k,v in x.items()}
  if isinstance(x,list):return [encode_detail(v) for v in x]
  return enc(x)
 return status,dict(constants={k:enc(v) for k,v in constants.items()},
  scalar_coefficient_derivative_bounds={n:[enc(v) for v in row] for n,row in sc.items()},
  connection_tensor_derivative_bounds=[enc(v) for v in dc],
  derivative_accounting=encode_detail(detail),cells=cells,cell_count=len(cells),
  exact_analytic_trial_additional_application_error=enc(arb(0)),
  implemented_caller_application_error=None,later_actual_minus_reference_application_error=None,signed_Gt=None,
  scope=dict(actual_seed_fractional_comparison=True,exact_analytic_trial_reference_residual_enclosed=True,
   actual_trial_H5_assumed=False,implemented_first_actual_inverse=False,complete_signed_quantum_correction=False,
   original_admission=False,nonlinear_existence=False,compatibility=False,physical_matching=False,
   validity_noise=False,all_sector_stability=False,construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());assert digest(Path(sys.executable))==rt['python_sha256']
 pkg=rt['packages']['flint'];assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
  python_sha256=rt['python_sha256'])
 try:
  rec['status'],rec['result']=execute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if rec['status']=='failed':raise SystemExit(1)
