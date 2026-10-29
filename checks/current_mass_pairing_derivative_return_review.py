"""Author proof review, original-norm ledger and complete spectral checks."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_mass_pairing_derivative_return as m
arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
ROOT=m.ROOT
digest=lambda p:sha256(p.read_bytes()).hexdigest()
path=ROOT/'evidence'/(m.STEM+'.json');raw=json.loads(path.read_bytes())
assert raw['status']=='massive_and_LLL_pairings_enclosed_full_mass_residual_open'
assert raw['input_sha256_before']==raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8']==(ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8']==(ROOT/'checks/current_mass_pairing_derivative_return.py').read_text(encoding='utf-8')
protected=json.loads(m.BASELINE.read_bytes())
assert all(digest(Path(p))==h for p,h in protected.items())
assert digest(Path(sys.executable))==m.runtime['python_sha256']
for pkg in m.runtime['packages'].values():
    assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
ctx.prec=768
r=raw['result'];c={k:dec(x) for k,x in r['constants'].items()}
sp=m.m.sp
t=sp.Symbol('t',positive=True);a=t*(1-t);P=sp.Integer(1)
for j in range(4):
    g=sp.exp(-1/a)
    assert sp.simplify(sp.diff(g,t,j)-g*P/a**(2*j))==0
    P=sp.expand(a*a*sp.diff(P,t)+sp.diff(a,t)*(1-2*j*a)*P)
# Geodesic chain rule, with z_ddot=Q'/2, verifies all physical coefficients.
z=sp.Symbol('z');u=sp.Symbol('zdot');Q=sp.Function('Q')(z);hh=sp.Function('h')(z)
D=lambda expr:sp.diff(expr,z)*u+sp.diff(expr,u)*sp.diff(Q,z)/2
expr=hh
jets=[]
for j in range(1,5):
    expr=sp.expand(D(expr));jets.append(expr)
fourth=(u**4*sp.diff(hh,z,4)+3*u*u*sp.diff(Q,z)*sp.diff(hh,z,3)
    +(sp.Rational(3,4)*sp.diff(Q,z)**2+2*u*u*sp.diff(Q,z,2))*sp.diff(hh,z,2)
    +(sp.diff(Q,z)*sp.diff(Q,z,2)/4+u*u*sp.diff(Q,z,3)/2)*sp.diff(hh,z))
assert sp.expand(jets[3]-fourth)==0
assert r['whole_cell_count']==1536 and r['cutoff_cell_count']==1024
assert c['cutoff_normalization']>0
assert c['gain_two_split']<arb(1)/400

q=arb('.338')
previous=m.m.constants('current-mass-pairing-integrable-return')
adj=m.m.constants('current-filtered-mass-adjoint-bound-attempt2')
bridge=m.m.constants('current-filtered-test-mass-trace-attempt2')
green=m.m.constants('current-filtered-mixed-green-residual')
alpha=previous['quantum_action_alpha_q'];Kg=adj['actual_filtered_mass_functional_norm_upper']
H1,H2,Bf=[c[k] for k in ['mass_shape_first_flow_L2_upper','mass_shape_second_flow_L2_upper','mass_shape_fractional_13_over_4_upper']]
N10,N11,N20,N21=[c[k] for k in ['actual_trial0_first_flow_norm_upper','actual_trial1_first_flow_norm_upper','actual_trial0_second_flow_norm_upper','actual_trial1_second_flow_norm_upper']]
eps=c['gain_two_epsilon'];kg=c['complete_gain_two_gradient_coefficient']
G=lambda a,b:(eps*b*b/3+kg*a*a).sqrt()
gain=alpha*Kg*G(H1,H2)*(q*G(N10,N20)+G(N11,N21))
principal=2*alpha*Kg*(q*N10+N11)*(H1+Bf)
assert gain.overlaps(c['complete_gain_two_mass_block_upper'])
assert principal.overlaps(c['complete_principal_mass_block_upper'])
total=principal+gain+c['complete_LLL_mass_block_upper']+previous['complete_zero_order_mass_block_upper']
assert total.overlaps(c['complete_massive_and_LLL_mass_pairing_upper'])
conv=bridge['full_E_to_retained_Et_norm_comparison'];eb=bridge['full_to_retained_j_energy_error']
remaining=(arb('1.49e34')-eb)/conv-total
assert remaining.overlaps(c['conditional_full_E_account_after_bridge_and_shape'])
assert (green['pref_t']*conv*total).overlaps(c['shape_block_original_bolt_charge'])
assert r['numerical_contract_pass']==bool(total<arb('1e34'))
assert r['complete_conditional_account_has_positive_remainder']==bool(remaining>0)
# Compare with the earlier first-derivative coefficient using the SAME
# fields and unchanged majorant X1. This is a sensitivity of estimators,
# not a lower bound on the actual source or residual.
oldc=previous['actual_order_one_reconstruction_constant']
oldN10=(2*arb.pi()).sqrt()*oldc*previous['trial0_order_one_auxiliary_norm_upper']
oldN11=(2*arb.pi()).sqrt()*oldc*previous['trial1_order_one_auxiliary_norm_upper']
oldgain=alpha*Kg*G(H1,H2)*(q*G(oldN10,N20)+G(oldN11,N21))
for name in ['actual_local_and_scalar_mass_residual','actual_scalar_return_bolt_value',
             'signed_r0_mass_pairing','signed_r1_mass_pairing','actual_signed_mass_vector_norm',
             'actual_mass_free_residual','complete_application_error','actual_Dj','actual_signed_Gt']:
    assert r[name] is None
for name in ['all_local_Weyl_scalar_and_classical_rows_included','complete_signed_mass_term_certified',
             'complete_signed_response_certified','complete_bolt_allowance','available_physical_margin',
             'original_P32_defect','actual_kappa_m','complete_original_admission','nonlinear_neighborhood_transfer',
             'nonlinear_existence','exact_compatibility','physical_matching','semiclassical_validity_noise',
             'all_sector_stability','construction_obstruction','full_goal_complete']:
    assert r['scope'][name] is False

ordinary=[
    'The exact f0/f1 trial, actual physical paired mass/log domain, original clock/temperature/spin/scale and fixed common regulator remain. Auxiliary spectral inequalities define and bound the original trial rather than selecting physical matter modes.',
    'The same specified normalized-bump mass transition is used. Its integral and derivatives are bounded on entire cells, including both flat endpoints. Cumulative bounds enclose the actual cutoff value. The original graph cutoff and mass conversion cost are not changed.',
    'The original primitive J_m is enclosed by the .15/.17 Q comparison, and all other jets follow its exact differential equation. Multiplying complete signed jets precedes absolute values, so derivative cutoff terms are present.',
    'The physical geodesic formulas follow exactly from z_ddot=Qprime/2. The angular cross term for ||P2h|| squared is 2pi A B, and its absolute value is retained. The H1 integral is the original pi int Q hprime squared, with no normalization improvement.',
    'The paired cubic is C2 with a third-derivative jump. The accepted seam-crossing difference argument and one-sided fourth jets give the numerical 13/4 bound. H2 core approximation and the complementary principal norms identify the pairings with the common-regulator form extension.',
    'The rational fj bound uses muY_nom<twice_a_nom and the full factor 1+2c_nom/muS_nom; it is not rounded down to two. The actual gamma and nominal source gamma remain distinct. Plancherel includes the entire auxiliary spectrum and actual averaging in the second derivative.',
    'The Fref bounds keep its complete inverse denominator through the proved Neumann ratio. Integrating powers of 1+p gives the exact constants in X11 and X21; a high-frequency asymptotic is not substituted.',
    'At derivative order one, only A grad f needs the source H1 norm. Every other product uses the complete auxiliary resolvent in L2, its gradient bound 1/(2sqrt(muY)) and Hessian bound. The exact D=3Qsecond-zQthird controls the full Volterra term. Both York rows and finite returns remain.',
    'For radial scalars the original SM first derivative is pi times the physical gradient integral. The actual Q/Q2<2 conversion gives sqrt(2pi) in (7). The accepted second-derivative comparison retains connection terms and applies at input order two, without an H5 assumption.',
    'The complete principal multiplier is bounded on complementary orders 13/4 and one. The actual mass fractional norm and source first derivative are now both quantitative.',
    'The gain-two form is polarized as a PRODUCT of combined norms. Its selected epsilon and logarithmic gradient cost, finite second/curvature contacts and all long geodesic returns are charged. No separated cross-product shortcut is used.',
    'The accepted zero-order block is reused unchanged. The fixed-domain LLL uses the same conformal-action bound. The local Weyl/anomaly and classical/scalar rows of equation (1) remain separate unknowns; this is the complete massive/LLL shape action, not the full coupled residual.',
    'Both mass residuals enter with q and one from ||(a0 T+a1 I)g||. The reference-adjoint error is already inside K_g and is not added or omitted separately.',
    'The full-to-retained response bridge and norm comparison remain in the original ledger. A positive conditional remainder is not an available physical margin or a bound for k_m, mass-free response or unknown numerical applications. All original nonlinear and physical obligations remain open.'
]
out=ROOT/'evidence'/(m.STEM+'-verification.json');assert not out.exists()
rec=dict(status='complete_massive_LLL_pairing_and_open_full_ledger_reviewed',
    reviewed_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(path),
    review_source_sha256=digest(Path(__file__)),frozen_inputs_checked=len(raw['input_sha256_after']),
    protected_paths_checked=len(protected),ordinary_review=ordinary,
    sensitivity=dict(previous_global_first_derivative_gain_estimator=enc(oldgain),
        previous_derivative_gain_estimator_over_complete_target=enc(oldgain/arb('1.49e34')),
        previous_over_attained_gain_estimator=enc(oldgain/gain)),
    numerical_contract_pass=r['numerical_contract_pass'],full_mass_residual_certified=False,
    original_admission=False,full_goal_complete=False,
    review_type='author ordinary proof and exact symbolic/interval verification',
    independent_specialist_review=False,lean_verification=False,
    review_source_utf8=Path(__file__).read_text(encoding='utf-8'))
out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in rec.items() if k not in ['ordinary_review','review_source_utf8']},indent=2))
