"""Author review and the next COMPLETE heat contract, not a trace theorem.

For the proposed alternative let h(t)[u,V eta] be the full time-reduced
heat Hessian minus its common a0/t and t*a2 terms, including their paired
boundary representative. If, on the actual fixed input and all unit Et
tests, |h(t)| <= Kheat*M_u*CV*t^(1+alpha), Theta(t)<=4C/t implies

  |c int_0^T Theta(t)h(t)dt/t|
      <=4cC*Kheat*M_u*CV*T^alpha/alpha, c=1/(B2 beta).

Take alpha=1/4, T=.001 and M_u=the accepted L2 upper plus the accepted
flow order-9/2 upper. This numerical scale bounds the lower flow norms;
it does NOT establish the required trace comparison. Assign HALF of the
conditional remaining residual account to this integral. The other half
must cover all local reallocations, positive-time terms, reference
comparison and actual application errors. None is assumed free or zero.
The resulting Kheat ceiling is a research target, not a proved constant.
This complete matched estimate would replace BOTH excessive cap contact
and gain-two reports together; it is not an isolated contact optimization.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-source-adapted-complete-response-test'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
def digest(p):return sha256(p.read_bytes()).hexdigest()
ctx.prec=512
receipt=ROOT/'evidence'/(STEM+'.json');raw=json.loads(receipt.read_bytes())
assert raw['status']=='complete_fixed_trial_estimator_excessive_signed_evaluation_open'
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_source_adapted_complete_response_test.py').read_text(encoding='utf-8')
protected=json.loads(BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
assert digest(Path(sys.executable))==rt['python_sha256']
pkg=rt['packages']['flint'];assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
res=raw['result'];v={k:dec(x) for k,x in res['constants'].items()}
d={k:dec(x) for k,x in res['debits'].items()}
assert set(d)=={'local_gain_two','finite_second_contact','global_row_subtracted_gluing',
 'actual_constant_contact_cap','actual_constant_contact_gluing','known_curvature_polynomial',
 'shared_vacuum_uncertainty','retained_energy_comparison'}
assert sum(d.values(),arb(0)).overlaps(v['complete_actual_minus_reference_upper'])
assert v['real_Born_ratio']<arb('.161') and v['complex_Born_ratio']<arb('.378')
assert v['real_squared_row_prefactor']<8 and v['complex_squared_row_prefactor']<32
# Independent cap L1 guard includes covariant endpoint/volume factors.
cap_L1=(arb('1.026')**3/arb('.974')).sqrt()*4/(1-4*arb('.026'))
assert cap_L1<10
# The constant-probe factors are rank 2, actual inverse 4/(3M),
# and Cauchy's 4M derivative plus M tadpole.
independent_pglue=arb(2)*4/3*5*32*260**2*4*792**2
assert independent_pglue.overlaps(v['constant_contact_gluing_coefficient'])
assert v['complete_actual_minus_reference_upper']>arb('1.66e44')
assert v['other_debits_without_cap_upper']>arb('8.53e39')
assert v['dominant_cap_fraction']>arb('.9999')
assert v['global_gluing_coefficient_improvement']>arb('1.02e11')
assert v['constant_contact_gluing_coefficient_improvement']>arb('4.56e10')
assert v['old_to_new_complete_upper_ratio']>121
assert res['actual_signed_nonzero_center'] is None and res['numerical_application_error'] is None
assert not any(res['scope'][k] for k in ['complete_actual_first_inverse','original_admission',
 'nonlinear_existence','compatibility','physical_matching','validity_noise','all_sector_stability',
 'construction_obstruction','goal_complete'])
praw=json.loads((ROOT/'evidence/current-global-green-trial-enclosure.json').read_bytes())
p={k:dec(x) for k,x in praw.get('runs',[praw])[-1]['result']['trial_parameters'].items()}
frac=json.loads((ROOT/'evidence/current-source-adapted-fractional-return.json').read_bytes())['result']['constants']
Mu=v['source_L2_upper']+dec(frac['actual_fractional_flow_norm_upper'])
alpha=arb(1)/4;T=arb('.001');account=arb('4.308597e32');allow=account/2
response_per_K=4*p['C']*Mu*v['test_weighted_L2_constant']*T**alpha/(alpha*p['B2']*p['beta'])
Kceiling=allow/response_per_K
entry_ratios={k:enc(x/account) for k,x in d.items()}
out=ROOT/'evidence'/(STEM+'-verification.json');assert not out.exists()
review=dict(status='complete_estimator_reviewed_excessive_no_signed_center',
 reviewed_utc=datetime.now(timezone.utc).isoformat(),review_type='author ordinary-argument and arithmetic/provenance review',
 independent_specialist_review=False,lean_verification=False,
 receipt_sha256=digest(receipt),review_source_sha256=digest(Path(__file__)),
 frozen_inputs_checked=len(raw['input_sha256_after']),protected_paths_checked=len(protected),
 cap_L1_independent_guard=enc(cap_L1),entry_over_CONDITIONAL_remaining_account=entry_ratios,
 proposed_complete_heat_contract=dict(alpha='1/4',T='.001',
  source_flow_scale_upper=enc(Mu),response_error_per_Kheat=enc(response_per_K),
  Kheat_ceiling_for_HALF_conditional_account=enc(Kceiling),
  attained_Kheat=None,actual_short_time_majorant_proved=False,
  remaining_half_account_certified=False,actual_application_errors_certified=False),
 ordinary_review=['The new resolvent rows retain both endpoint buffers and the actual global inverse.',
  'All Born insertions and the full angular count are summed; actual thermal returns stay inside A.',
  'The complete constant probe includes its tadpole, common finite restoration, and row correction.',
  'The accepted bounded-form extension, not finite-mass analyticity, defines the fixed-input response.',
  'Eliminating the cap contact alone leaves an excessive complete remaining report.',
  'The heat comparison is a prospective sufficient contract, not an established theorem.'])
out.write_bytes((json.dumps(review,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=review['status'],Kheat_ceiling=enc(Kceiling)['ball'],
 per_K=enc(response_per_K)['ball'],entry_ratios={k:x['ball'] for k,x in entry_ratios.items()},
 inputs=review['frozen_inputs_checked'],protected=review['protected_paths_checked']),indent=2))
