"""Author review of the structured integrability gain and its limited budget scope."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_mass_pairing_integrable_return as m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
ROOT=m.ROOT
digest=lambda p:sha256(p.read_bytes()).hexdigest()
path=ROOT/'evidence'/(m.STEM+'.json')
raw=json.loads(path.read_bytes())
assert raw['status']=='structured_L1_return_and_zero_order_mass_block_bounded_full_response_open'
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_mass_pairing_integrable_return.py').read_text(encoding='utf-8')
protected=json.loads(m.BASELINE.read_bytes())
assert all(digest(Path(p))==h for p,h in protected.items())
assert digest(Path(sys.executable))==m.runtime['python_sha256']
for pkg in m.runtime['packages'].values():
    assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())

ctx.prec=768
r=raw['result'];c={k:dec(v) for k,v in r['constants'].items()}
v={k:m.rational(x) for k,x in m.m.result('current-correlated-signed-trial-test-attempt3')['exact_nominal_parameters'].items()}
# Exact rational factorization on the entire auxiliary spectrum, not a
# finite approximation to the physical field operator or quantum tower.
sp=m.m.sp
p=sp.Symbol('p',nonnegative=True)
ss=lambda x:sp.Rational(x.numerator,x.denominator)
vs={k:ss(x) for k,x in v.items()}
fj=(p+vs['muY'])*(p+vs['muS']+2*vs['scalar_c'])/((p+vs['twice_a'])*(p+vs['muS'])*(p+vs['muT']))
factor=(1+(vs['muY']-vs['twice_a'])/(p+vs['twice_a']))*(1+2*vs['scalar_c']/(p+vs['muS']))/(p+vs['muT'])
assert sp.cancel(fj-factor)==0
h=vs['d1']/2
Delta=h*h-vs['d0']
assert sp.expand((p+h)**2-Delta-(p*p+vs['d1']*p+vs['d0']))==0
assert sp.expand((vs['Cv']-vs['delta']*h)+vs['delta']*(p+h)-(vs['Cv']+vs['delta']*p))==0

adj=m.m.constants('current-filtered-mass-adjoint-bound-attempt2')
transfer=m.m.constants('current-filtered-test-mass-trace-attempt2')
q=arb('.338');kg=adj['actual_filtered_mass_functional_norm_upper']
combination=kg*(q*c['complete_trial0_shape_L1_upper']+c['complete_trial1_shape_L1_upper'])
glue=combination*c['quantum_action_alpha_q']*arb('6e24')*c['actual_mass_shape_Linfinity_upper']
contact=combination*2*c['constant_contact_Lq']*arb('1e20')*c['actual_weighted_contact_factor']
assert glue.overlaps(c['complete_gluing_block_upper'])
assert contact.overlaps(c['complete_constant_contact_block_upper'])
assert (glue+contact).overlaps(c['complete_zero_order_mass_block_upper'])
assert c['complete_zero_order_mass_block_upper']<arb('1e28')
assert c['zero_order_block_fraction_of_conditional_target']<arb('1e-6')
assert c['rejected_old_estimator_lower_over_new_zero_order_upper']>arb('1e12')

# Independently assemble the first derivative comparison from all terms
# of the accepted actual-Q reconstruction, including its Volterra term.
fr=m.m.result('current-source-adapted-fractional-return')
sc={k:[dec(x) for x in row] for k,row in fr['scalar_coefficient_derivative_bounds'].items()}
dc=[dec(x) for x in fr['connection_tensor_derivative_bounds']]
params,geo,bs=m.fr.load()
_,_,detail=m.fr.full_map_constants(params,sc,dc)
S1=detail['reconstructed_reference_derivative_constants'][1]
# full_map_constants takes upper ENDPOINTS at each product. Reassembling
# from decimal-enclosed input constants gives a different upper estimator,
# so overlap with the old thin endpoint is not a valid requirement. The
# original accepted S1 remains the bound used by the ordinary proof;
# reconstructing from its original receipt must reproduce that usage.
accepted_S1=dec(fr['derivative_accounting']['reconstructed_reference_derivative_constants'][1])
assert (arb(10).sqrt()*accepted_S1).overlaps(c['actual_order_one_reconstruction_constant'])
reassembly_relative_difference=abs(S1-accepted_S1)/accepted_S1
assert reassembly_relative_difference<arb('1e-55')
assert S1<934 and accepted_S1<934
conv=transfer['full_E_to_retained_Et_norm_comparison']
bridge=transfer['full_to_retained_j_energy_error']
remaining=(arb('1.49e34')-bridge)/conv-glue-contact
assert remaining.overlaps(c['conditional_full_E_account_after_bridge_and_zero_order'])
assert (remaining/c['principal_pairing_coefficient_Ageo']).overlaps(c['conditional_mass_shape_ceiling_ALL_other_errors_zero'])
assert remaining>0
for key in ['mass_shape_fractional_norm','principal_mass_block','complete_derivative_and_scalar_mass_remainder',
            'signed_r0_mass_pairing','signed_r1_mass_pairing','actual_signed_mass_vector_norm',
            'mass_free_residual','complete_application_error','actual_Dj','actual_signed_Gt']:
    assert r[key] is None
for key in ['full_signed_mass_term_certified','signed_response_certified','complete_bolt_allowance',
            'available_physical_margin','original_P32_defect','actual_kappa_m','complete_original_admission',
            'nonlinear_neighborhood_transfer','nonlinear_existence','exact_compatibility','physical_matching',
            'semiclassical_validity_noise','all_sector_stability','construction_obstruction','full_goal_complete']:
    assert r['scope'][key] is False

ordinary=[
    'The physical paired domain, mass/log phases, exact ultraviolet logarithm, common regulator, complete angular/thermal sums, all contacts, cutoff transitions and regular bolt remain. Auxiliary H defines the unchanged trial; it is not a replacement physical operator.',
    'The original mass cutoff endpoints .1/.2 were specified but no intervening profile was fixed. The normalized smooth bump now fixes it explicitly. The zero-order estimates hold for every such 0<=chi<=1; neither improvement nor original conversion cost is attributed to choosing it.',
    'Positive Dirichlet exhaustion on the auxiliary disk yields the original Friedrichs resolvent and its L1 contraction 1/mu. L1 and L2 realizations agree on the trial. The rational identities and full Neumann series bound apply to the complete auxiliary spectrum.',
    'Actual muY is retained in reconstruction and is not canceled against a nominal pole. Original dz=gamma dx, and the nominal gamma inside A0, are both retained.',
    'For r in L1 with Hr in L1, the flux has a finite limit. Its nonzero limit would give a nonintegrable 1/x tail, so it vanishes. Integrating only over the finite physical length bounds both the flux and Volterra terms without an unproved first moment.',
    'Both York rows give the whole reconstruction A f+B r+C Hflux+E W. Its triangle estimate includes every term and both finite endpoints; it bounds the actual trial shape in L1.',
    'The accepted gluing proof bounds BOTH absolute kernel rows by 3e24 C. Row subtraction adds the diagonal action and gives the L-infinity operator bound 6e24 C. Pairing against the actual L1 trial is therefore justified, unlike inferring an L-infinity bound from an L2 norm.',
    'The mass cubic compensates the constant-contact 1/Y factor on its complete support, including the transition. The constant contact and row-subtracted gluing are the two separate blocks of the existing common scheme; no row term is dropped or counted twice.',
    'These fixed-input kernel pairings have integrable majorants and agree with the original form extension. The zero-order argument does not claim a new nonlinear metric-variation theorem or control the ultraviolet derivative block.',
    'The prescribed averaged source has a finite first auxiliary derivative norm. The complete spectral integral is bounded using the density and the actual bump Plancherel norm. No physical quantum mode is removed. The point source cannot replace this averaged input at order one.',
    'At input order one the complete reconstruction uses at most two derivatives of R_muY f. Its residual spectral multiplier is bounded by one; the accepted S1 product coefficient and sqrt(10) metric conversion therefore apply without requiring an H5 input.',
    'The explicit mass lift is C2 across the actual seam with bounded one-sided fourth jets. The accepted seam-crossing difference proof gives Dom(|P|^(13/4)), sufficient for pairing with the order-one source and the complete 17/4 multiplier. No useful numerical mass-shape norm is supplied.',
    'In York variables f_m=b h_m-a s_m and k_m=a^-1 f_m. The local scalar row is supported with Phi_m; its bolt average is a separate unknown. The two exact residual identities preserve that inverse and both canonical h pairings. An earlier different scalar source cannot be substituted.',
    'The full-to-retained norm conversion multiplies the entire new full-E block. The remaining account and M_h ceiling grant every other unknown debit zero and are conditional only. There is no complete signed mass, Dj, Gt, bolt, admission or stability certificate.'
]
out=ROOT/'evidence'/(m.STEM+'-verification.json')
assert not out.exists(),'Preserve immutable receipts'
failed_path=ROOT/'evidence'/(m.STEM+'-verification-failed.json')
failed=json.loads(failed_path.read_bytes())
assert failed['status']=='failed'
assert failed['arithmetic_receipt_sha256']==digest(path)
assert failed['ordinary_argument_sha256']==digest(ROOT/'notes'/(m.STEM+'.md'))
rec=dict(status='structured_L1_gain_and_limited_budget_contribution_reviewed',
    reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(path),
    review_source_sha256=digest(Path(__file__)),frozen_inputs_checked=len(raw['input_sha256_after']),
    protected_paths_checked=len(protected),ordinary_review=ordinary,
    exact_whole_spectrum_factorizations=True,complete_first_derivative_product_reassembled=True,
    derivative_reassembly_relative_difference=enc(reassembly_relative_difference),
    derivative_bound_basis='Unchanged accepted S1 from its frozen original enclosure. Reassembly from decimal-enclosed coefficients is a nearby upper estimate, not the identical thin endpoint.',
    preserved_failed_review=dict(path=str(failed_path.relative_to(ROOT)),sha256=digest(failed_path),
        cause='Invalid overlap requirement between upper endpoints produced before and after decimal input enclosure. Original derivative bound and ordinary proof unchanged.'),
    complete_zero_order_mass_block_bounded=True,complete_signed_mass_term_certified=False,
    available_physical_margin=False,full_goal_complete=False,
    review_type='author ordinary argument and exact symbolic/interval arithmetic',
    independent_specialist_review=False,lean_verification=False,
    review_source_utf8=Path(__file__).read_text(encoding='utf-8'))
out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in rec.items() if k not in ['ordinary_review','review_source_utf8']},indent=2))
