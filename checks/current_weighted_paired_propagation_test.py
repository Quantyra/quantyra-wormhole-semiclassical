"""Paired-domain weight and Hspin propagation; quantum trace bound unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-weighted-paired-propagation-test'
INPUTS=['notes/'+STEM+'.md','checks/current_weighted_paired_propagation_test.py',
 'checks/current_quartic_ritz_rank_test.py','notes/current-rough-paired-variation.md',
 'notes/stationary-massive-boundary-layer.md','notes/current-disk-heat-propagation.md',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-localized-short-time-trial.json',
 'evidence/current-localized-covariant-jet-test.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(name):
 d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
 r=d.get('runs',[d])[-1];assert r['status']!='failed'
 return r['result']
def execute():
 ctx.prec=512
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 h,L,gamma,beta,ep=[p[k] for k in ['h','L','gamma','beta','epsilon_seed']]
 delta=h/(1+h);B=gamma+L;w0=gamma/B
 assert 0<w0<1 and 2*beta-4*abs(ep)>0
 assert (6*beta+p['u']+ep).contains(1)
 ratio=delta/gamma
 alpha=3*arb(max(arb(1).upper(),ratio.upper()))
 t0=arb('.001');base=(alpha**2*t0).exp();time=2*(arb(3)/2*alpha**2*t0).exp()
 robin=3*beta.sqrt()/B
 Wnorm=(gamma*(1-w0**5)/5).sqrt()
 trial=result('current-localized-short-time-trial')['constants']
 Dabs=dec(trial['conditional_bulk_absolute_ceiling_OTHER_ZERO'])
 Wceiling=Dabs/Wnorm
 constants=dict(alpha_upper=alpha,original_boundary_relative_Robin_term=robin,
  W_boundary_value=w0**3,original_W_L2=Wnorm,
  crude_Hspin_similarity_bound=w0**-3,
  correct_Hspin_real_time_bound=base,correct_Hspin_Dt_bound=time,
  Hspin_bound_improvement_over_crude=w0**-3/base,
  conditional_complete_weighted_density_ceiling_OTHER_ZERO=Wceiling,
  norm_ratio_H8_over_W=dec(result('current-localized-covariant-jet-test')['constants']['actual_full_H8_L2'])/Wnorm)
 assert robin>1 and base<arb('1.01') and time<arb('2.03')
 assert all(v.is_finite() for v in constants.values())
 return dict(precision_bits=ctx.prec,t0='0.001',constants={k:enc(v) for k,v in constants.items()},
  ordinary_domain_test='W preserves Dom A but does not preserve the original Neumann part of Dom A squared',
  operator_space='original orbit L2 spinors, not metric E_t or a diagonal-density space',
  complete_weighted_kernel_constant=None,complete_quantum_remainder=None,complete_signed_first_residual=None,
  scope=dict(correct_form_similarity_defined=True,all_thermal_and_radial_modes=True,
   uniform_mass_parameter_estimate=True,both_weight_directions_bounded=True,
   weighted_diagonal_or_Schur_estimate=False,angular_quantum_trace_evaluated=False,
   complete_quantum_certificate=False,new_physical_admission_margin=False,
   nonlinear_existence_matching_validity_stability=False,construction_obstruction=False,goal_complete=False))
if __name__=='__main__':
 out=ROOT/'evidence'/(STEM+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
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
  rec.update(status='weighted_form_propagation_bounded_complete_quantum_trace_unresolved',
   input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
