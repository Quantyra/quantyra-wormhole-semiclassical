"""Actual-domain localized trial; complete quantum remainder remains unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc

ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-short-time-trial'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_localized_short_time_trial.py',
 'checks/current_quartic_ritz_rank_test.py',
 'notes/current-local-reference-range-test.md',
 'notes/current-smeared-heat-projection-domain.md',
 'notes/current-disk-heat-propagation.md',
 'notes/current-positive-time-quantum-operator.md',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-positive-time-quantum-operator.json',
 'evidence/current-first-quantum-solve-attempt.json',
 'evidence/current-first-inverse-source-screen.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p): return sha256(p.read_bytes()).hexdigest()
def result(name):
 d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
 r=d.get('runs',[d])[-1]
 assert r['status']!='failed'
 return r['result']
def mono(n,c=1,k=0): return {(n,k):arb(c)}
def clean(p): return {k:v for k,v in p.items() if not v.is_zero()}
def add(*ps):
 out={}
 for p in ps:
  for k,v in p.items(): out[k]=out.get(k,arb(0))+v
 return clean(out)
def scale(p,c): return clean({k:v*c for k,v in p.items()})
def mul(p,q):
 out={}
 for (n,k),v in p.items():
  for (m,j),u in q.items():
   key=(n+m,k+j);out[key]=out.get(key,arb(0))+v*u
 return clean(out)
def power(p,n):
 out=mono(0)
 for _ in range(n):out=mul(out,p)
 return out
def diff(p):
 out={}
 for (n,k),v in p.items():
  if n:out=add(out,mono(n-1,v*n,k))
  if k:out=add(out,mono(n-1,v*k,k-1))
 return out
def val(p,w):
 return sum((v*w**n*w.log()**k for (n,k),v in p.items()),arb(0))
def primitive_from(p,w0):
 out={}
 for (n,k),v in p.items():
  assert k==0
  out=add(out,mono(0,v,1) if n==-1 else mono(n+1,v/(n+1)))
 return add(out,mono(0,-val(out,w0)))
def moment(n,k,w0):
 if n==-1:return -w0.log()**(k+1)/(k+1)
 a=n+1;pw=w0**a
 out=(1-pw)/a
 for j in range(1,k+1):out=-pw*w0.log()**j/a-j*out/a
 return out
def l2(p,w0,gamma):
 square=mul(p,p)
 integral=gamma*sum((v*moment(n-2,k,w0) for (n,k),v in square.items()),arb(0))
 assert integral>0,integral
 return integral.sqrt()
def serialize(p):
 return [dict(power=n,log_power=k,coefficient=enc(v)) for (n,k),v in sorted(p.items())]

def execute():
 ctx.prec=512
 gp=result('current-global-green-trial-enclosure')['trial_parameters']
 p={k:dec(v) for k,v in gp.items()}
 L,gamma,a,beta,C,B2=[p[k] for k in ['L','gamma','a','beta','C','B2']]
 mp=result('current-global-mass-conversion')['parameters']
 u,ep,tau=[dec(mp[k]) for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+ep,-2*ep,ep]
 qm=add(mono(-1),mono(0,-1))
 qb=add(scale(qm,2),scale(power(qm,2),a))
 exact_coefficients={1:2*gamma,2:a}
 for j in range(3,7):
  sj=(-1)**j*sum((comb(i,j)*Q[i]*L**(i-j) for i in range(j,7)),arb(0))
  exact_coefficients[j]=sj
  qb=add(qb,scale(power(qm,j),sj*gamma**(j-2)))
 z=add(mono(0,L+gamma),mono(-1,-gamma));w0=gamma/(gamma+L)
 dw=lambda f:mul(mono(2),diff(f))
 delta=lambda f:dw(mul(qb,dw(f)))
 Qp=scale(dw(qb),gamma);Qpp=dw(dw(qb));R=scale(Qpp,-1)
 M=add(mul(z,Qpp),scale(Qp,-4))
 If=primitive_from(mul(M,mono(2)),w0)
 A=dw(mul(mul(power(z,5),qb),mono(4)))
 psi=add(A,scale(mul(power(z,4),If),-1))
 v=add(scale(A,-arb(3)/2),mul(power(z,4),If))
 un=add(scale(psi,arb(1)/2),scale(mul(mul(power(z,4),qb),mono(4)),gamma))
 vnorm=l2(v,w0,gamma)
 jets=[vnorm];cur=v
 for j in range(1,5):
  cur=delta(cur);jets.append(l2(cur,w0,gamma))
 H=lambda f:add(delta(f),mul(R,f))
 dv=delta(v);dR=delta(R);R2=mul(R,R)
 O6=add(scale(add(H(delta(H(v))),mul(dR,H(v)),delta(mul(v,dR))),-arb(1)/105),
  scale(add(scale(mul(mul(R2,R),v),16),scale(mul(R2,dv),24),
            scale(delta(mul(R2,v)),24),scale(delta(mul(R,dv)),24)),-arb(1)/1008))
 O6norm=l2(O6,w0,gamma)
 plus=result('current-positive-time-quantum-operator')['constants']
 CLt=dec(plus['CL_t']);k=dec(plus['total_operator_difference_upper'])
 rho=dec(plus['conditional_short_residual_ceiling'])
 Ghat=dec(result('current-first-quantum-solve-attempt')['constants']['proved_source_dual_upper'])
 Kfirst=dec(result('current-first-inverse-source-screen')['constants']['first_residual_to_original_bolt_multiplier'])
 energy_lower=vnorm/CLt;t0=arb('.001')
 pref=4*C*CLt**2/(B2*beta)
 tau6=pref*t0*O6norm/(4*arb.pi()*vnorm)
 tau_req=(1-k)*rho/Ghat
 remaining=tau_req-tau6
 assert remaining>0
 cb=pref*arb(4)/7*t0**(arb(7)/4);cd=pref*t0**2/2
 constants=dict(original_v_L2=vnorm,trial_Et_lower=energy_lower,
  original_O6_L2=O6norm,O6_over_v_L2=O6norm/vnorm,
  physical_measure_gamma=gamma,global_to_normalized_density_cost=1/vnorm,
  attained_leading_column_upper=tau6,conditional_relative_column_screen=tau_req,
  leading_upper_fraction_of_conditional_screen=tau6/tau_req,
  leading_original_first_bolt_upper=Kfirst*Ghat*tau6/(1-k),
  conditional_boundary_relative_ceiling_OTHER_ZERO=remaining/cb,
  conditional_bulk_relative_ceiling_OTHER_ZERO=remaining/cd,
  conditional_boundary_absolute_ceiling_OTHER_ZERO=remaining/cb*vnorm,
  conditional_bulk_absolute_ceiling_OTHER_ZERO=remaining/cd*vnorm,
  bulk_relative_budget_per_Delta4_ratio=remaining/cd/(jets[4]/vnorm),
  boundary_budget_coefficient=cb,bulk_budget_coefficient=cd,
  seed_energy_shape_comparison=CLt)
 for j,n in enumerate(jets):
  constants['Delta'+str(j)+'_L2']=n
  constants['Delta'+str(j)+'_over_v_L2']=n/vnorm
 assert all(x.is_finite() for x in constants.values())
 return dict(precision_bits=ctx.prec,constants={k:enc(v) for k,v in constants.items()},
  exact_trial=dict(coordinate='w=gamma/(gamma+L-z)',w0=enc(w0),
   Q_over_gamma_squared=serialize(qb),psi=serialize(psi),Upsi=serialize(un),Vpsi=serialize(v)),
  contracted_jet_order=8,complete_quantum_remainder_Bv=None,complete_quantum_remainder_Dv=None,
  actual_Ritz_coefficient=None,complete_first_residual=None,signed_Gt=None,
  scope=dict(actual_global_trial_domain_proved=True,both_York_rows_retained=True,
   full_finite_Q_and_original_measure=True,whole_interval_trial_regularities_enclosed=True,
   complete_quantum_action_evaluated=False,new_complete_admission_margin=False,
   physical_norm_comparison_for_actual_solution=False,nonlinear_or_independent_phase_bound=False,
   existence_matching_validity_stability_proved=False,construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 path=ROOT/'evidence'/(STEM+'.json')
 assert not path.exists(),'Preserve immutable execution receipt.'
 before={p:digest(ROOT/p) for p in INPUTS}
 protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(status='trial_regularity_enclosed_complete_short_time_response_unresolved',
   input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 path.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 keys=['original_v_L2','original_O6_L2','Delta4_over_v_L2',
  'leading_upper_fraction_of_conditional_screen','conditional_bulk_absolute_ceiling_OTHER_ZERO',
  'conditional_bulk_relative_ceiling_OTHER_ZERO','bulk_relative_budget_per_Delta4_ratio']
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:rec['result']['constants'][k]['ball'] for k in keys if 'result' in rec}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
