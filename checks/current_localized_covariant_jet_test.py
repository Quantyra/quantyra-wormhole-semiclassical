"""Exact tensor recurrence on the actual trial; no quantum heat certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_localized_short_time_trial import arb,ctx,dec,enc,mono,add,mul,scale,diff,power,moment
ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-covariant-jet-test'
INPUTS=['notes/'+STEM+'.md','checks/current_localized_covariant_jet_test.py',
 'notes/current-localized-short-time-trial.md','checks/current_localized_short_time_trial.py',
 'checks/current_quartic_ritz_rank_test.py','evidence/current-localized-short-time-trial.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
def digest(p):return sha256(p.read_bytes()).hexdigest()
def read(items):return {(x['power'],x['log_power']):dec(x['coefficient']) for x in items}
def execute():
 ctx.prec=512
 old=json.loads((ROOT/INPUTS[5]).read_bytes())
 assert old['status']=='trial_regularity_enclosed_complete_short_time_response_unresolved'
 assert all(digest(ROOT/k)==v for k,v in old['input_sha256_after'].items())
 r=old['result'];q=read(r['exact_trial']['Q_over_gamma_squared']);v=read(r['exact_trial']['Vpsi'])
 w0=dec(r['exact_trial']['w0']);gamma=dec(r['constants']['physical_measure_gamma'])
 D=lambda f:mul(mono(2),diff(f));h=scale(D(q),arb(1)/2)
 tensors={():v};norms=[];sizes=[]
 qpowers=[power(q,i) for i in range(5)]
 for m in range(9):
  density={}
  for mask in range(2**m):
   if mask.bit_count()%2:continue
   component={}
   for pairs,coeff in tensors.items():
    paired={i for pair in pairs for i in pair}
    if any((mask>>i)&1 for i in range(m) if i not in paired):continue
    if any(((mask>>i)&1)!=((mask>>j)&1) for i,j in pairs):continue
    component=add(component,mul(coeff,qpowers[(m-2*len(pairs))//2]))
   density=add(density,mul(component,component))
  if m%2:density=mul(q,density)
  integral=gamma*sum((c*moment(n-2,k,w0) for (n,k),c in density.items()),arb(0))
  assert integral>0,(m,integral)
  norms.append(integral.sqrt());sizes.append(len(tensors))
  if m==8:break
  new={}
  for pairs,coeff in tensors.items():
   new[pairs]=add(new.get(pairs,{}),D(coeff))
   paired={i for pair in pairs for i in pair}
   for i in range(m):
    if i in paired:continue
    key=tuple(sorted(pairs+((i,m),)))
    new[key]=add(new.get(key,{}),mul(h,coeff))
  tensors=new
 assert norms[0].overlaps(dec(r['constants']['original_v_L2']))
 H8=sum(norms,arb(0));Hrel=H8/norms[0]
 bulk=dec(r['constants']['conditional_bulk_relative_ceiling_OTHER_ZERO'])
 constants=dict(actual_full_H8_L2=H8,actual_full_H8_over_v_L2=Hrel,
  conditional_complete_bulk_H8_constant_ceiling_OTHER_ZERO=bulk/Hrel,
  global_pointwise_to_normalized_density_cost=1/norms[0])
 for i,x in enumerate(norms):constants['covariant_'+str(i)+'_L2']=x
 return dict(precision_bits=ctx.prec,tensor_pairing_counts=sizes,
  constants={k:enc(v) for k,v in constants.items()},
  quantum_remainder_constant=None,complete_signed_first_residual=None,
  scope=dict(full_actual_trial_covariant_regularities_enclosed=True,
   arbitrary_energy_to_H8_comparison=False,full_physical_admission=False,
   complete_quantum_heat_remainder_certified=False,goal_complete=False))
if __name__=='__main__':
 out=ROOT/'evidence'/(STEM+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS}
 protected=json.loads(BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert before==after
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(status='full_trial_tensor_regularity_enclosed_heat_remainder_unknown',
   input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 keys=['actual_full_H8_L2','actual_full_H8_over_v_L2',
  'conditional_complete_bulk_H8_constant_ceiling_OTHER_ZERO']
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:rec['result']['constants'][k]['ball'] for k in keys if 'result' in rec}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
