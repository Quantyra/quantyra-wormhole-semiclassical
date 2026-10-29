"""Author review of the actual-domain endpoint/interpolation comparison."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_source_adapted_fractional_return as f
arb,ctx,dec,enc,J=f.arb,f.ctx,f.dec,f.enc,f.J
ROOT=f.ROOT;STEM=f.STEM
def digest(p):return sha256(p.read_bytes()).hexdigest()

ctx.prec=512
out=ROOT/'evidence'/(STEM+'-verification.json');assert not out.exists()
raw=ROOT/'evidence'/(STEM+'.json');d=json.loads(raw.read_bytes())
assert d['status']=='fractional_map_enclosed_exact_trial_screen_passes_numerical_applications_open'
assert d['input_sha256_before']==d['input_sha256_after']
for p,h in d['input_sha256_before'].items():assert digest(ROOT/p)==h,p
assert d['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
assert d['arithmetic_source_utf8']==(ROOT/'checks/current_source_adapted_fractional_return.py').read_text(encoding='utf-8')
assert digest(Path(sys.executable))==d['python_sha256']
baseline=json.loads(f.lr.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in baseline.items())
r=d['result'];c={k:dec(v) for k,v in r['constants'].items()}
assert c['C4']<arb('5.528986e8') and c['C5']<arb('5.219650e10')
assert c['actual_reconstruction_fractional_comparison']<arb('5.372092e9')
assert c['actual_fractional_flow_norm_upper']<arb('4.820204e58')
assert c['complete_reference_nonlocal_response']<arb('7.494891e26')
assert c['quantum_debit_over_remaining_account']<arb('1.739517e-6')
assert c['complete_exact_analytic_first_reference_residual']<arb('5.645453e32')
assert c['complete_exact_analytic_first_reference_residual']<arb('9.955e32')
assert c['reference_inverse_energy_error']<arb('5.658e32')
assert c['remaining_for_UNEXECUTED_residual_application_errors']>arb('4.308597e32')
assert c['known_coarse_actual_minus_reference_upper']>arb('2.02e46')
assert c['coarse_actual_difference_upper_over_complete_first_screen']>arb('2.02e13')
assert r['implemented_caller_application_error'] is None and r['signed_Gt'] is None
assert r['later_actual_minus_reference_application_error'] is None
assert dec(r['exact_analytic_trial_additional_application_error']).is_zero()
assert not r['scope']['actual_trial_H5_assumed']

p,old,bs=f.load();cells=r['cells']
assert dec(cells[0]['s0']).is_zero() and dec(cells[-1]['s1']).overlaps(p['L'])
for a,b in zip(cells,cells[1:]):assert a['s1']==b['s0']
sc={n:[dec(v) for v in row] for n,row in r['scalar_coefficient_derivative_bounds'].items()}
dc=[dec(v) for v in r['connection_tensor_derivative_bounds']]
for row in cells:
 for n in sc:
  for upper,cell in zip(sc[n],row['scalar_bounds'][n]):assert dec(cell).upper()<=upper.upper()
 for upper,cell in zip(dc,row['connection_bounds']):assert dec(cell).upper()<=upper.upper()

# Expanded connection-transfer polynomials independently check recurrence (5).
S=[dec(v) for v in r['derivative_accounting']['reconstructed_reference_derivative_constants']]
D0,D1,D2,D3=dc
T4=S[4]+6*D0*S[3]+(4*D1+11*D0**2)*S[2]+(D2+7*D0*D1+6*D0**3)*S[1]
T5=S[5]+10*D0*S[4]+(10*D1+35*D0**2)*S[3]+(5*D2+45*D0*D1+50*D0**3)*S[2]
T5+=(D3+11*D0*D2+7*D1**2+46*D0**2*D1+24*D0**4)*S[1]
assert c['C4'].overlaps(100*T4)
assert c['C5'].overlaps(arb(10)**(arb(5)/2)*T5)
assert c['actual_reconstruction_fractional_comparison'].overlaps((c['C4']*c['C5']).sqrt())

# Independently recover the connection from metric compatibility and torsion.
fixtures=[]
for value in [arb('.3'),arb(2),arb('1e20')]:
 x=J([value,arb(1)]+[arb(0)]*9)
 z,q2,R,q=f.geometry_jets(p,bs,x)
 radial=f.sqrtjet(q2);Rprime=radial*R.d();A2=q2.d()/(2*radial)
 b=Rprime/(2*R);cc=(1-R**2)*A2-R*Rprime/2
 checks=[radial*(1/R).d()+2*b/R,
         Rprime-2*R*b,
         A2*(1/R-R)-R*b-cc/R]
 assert all(t.c[0].contains(0) for t in checks)
 # Hessian x=(1+a*x)g2 and all smooth-core bolt identities are ordinary formulas.
 assert (q2.d()/2-(1+p['a']*x)).c[0].contains(0)
 fixtures.append(dict(x=enc(value),metric_compatibility=True,torsion_free=True,reference_Hessian_identity=True))

# Whole-reference residual and its coercivity conversion, without dropping unknowns.
first=f.lr.constants('current-first-quantum-solve-attempt')
e0=f.up(first['complete_geodesic_form_error'])
assert c['reference_inverse_energy_error'].overlaps(c['complete_exact_analytic_first_reference_residual']/(1-e0))
assert c['reference_inverse_energy_error']<first['conditional_common_stage_output_error_ceiling']
assert c['remaining_for_UNEXECUTED_caller_Et_errors']>arb('4.32e32')
trial_upper=(f.up(first['proved_source_dual_upper'])+c['complete_exact_analytic_first_reference_residual'])/(1-e0)
assert c['actual_trial_Et_upper'].overlaps(trial_upper)
assert c['known_coarse_actual_minus_reference_upper'].overlaps(e0*trial_upper)
reviews=dict(
 spaces='Auxiliary H Sobolev endpoint inputs define a linear actual finite-Q York reconstruction. Physical output is the original dz/averaged-direction norm of the complete paired flow generator P. No physical-mode truncation or normalization change.',
 endpoint_domain='Smooth auxiliary core inputs and their regular resolvents reconstruct smooth bolt fields. Their outer z^4 vanishing yields matching fourth flow jets on the C2,1 paired seam; the fifth weak flow derivative has no delta. Closure, not a smooth surrogate theorem, defines the endpoint maps.',
 model_derivatives='On curvature -a, parallel curvature gives commutator norm <=a sum(1+6j)<=4k^2 a. Scalar invariant differential contractions are polynomials in H; integration by parts and spectral calculus give the stated explicit derivative constants through order six.',
 actual_geometry='The metric ratio remains .1<R<2 and volume is unchanged. The exact connection components and all three reference covariant derivatives are bounded on whole cells; h^3 factoring retains the regular-bolt cancellation. Coordinate-free tensor norms and 10^(n/2) restore actual units.',
 field='All four terms of the exact two-row shape reconstruction are retained, including the whole Volterra return. Resolvent gain supplies the one extra derivative of r in the gradient term. Multinomial and Bell bounds apply to full covariant derivative tensors, not sampled radial derivatives.',
 interpolation='The same linear map is bounded at orders four and five. Spectral three-lines interpolation gives order 9/2, with Cflow=sqrt(C4 C5). It never asserts that the actual rho^4 log input is H5; the accepted H^(9/2) model-profile norm is used only after the actual output comparison is proved.',
 full_quantum_scope='The accepted Kgeo then encloses the complete reference nonlocal column, exact UV logarithm and all angular/thermal/flow returns, with common local contacts already charged. No actual sign is claimed. Actual-minus-reference quantum applications are not covered by this new number.',
 application_scope='All interval coefficient and bound-arithmetic errors are included for the exact analytic trial. Its additional numerical application is exactly absent by definition. Errors of a numerical realization/caller and later ordered operator applications remain UNKNOWN and nonzero unless separately certified.',
 norm_and_budget='The complete reference Et* residual is converted to an Et inverse error with 1/(1-e0). Original bolt coefficient products additionally retain that factor and the known repeated-source amplification. They are components, not full physical-output or admission bounds.',
 coarse_actual_estimator='The known general-form estimate for (I-A0)x_hat is >2.02e46 as an upper expression, >2.02e13 times the complete first-residual screen. This is an excessive fallback report, not a lower bound on the actual response or an obstruction. Its next signed evaluation must distinguish amplitude from numerical error.',
 decision='Retain the reconstruction/interpolation method and stop refining this now affordable quantum reference term. Next enclose the actual-minus-reference quantum action on the certified analytic trial, with full common-regulator source hypotheses and application errors.',
 limits='No complete actual inverse, Gt, original P32 defect, mass column, graph/physical output, nonlinear admission/existence, exact compatibility, physical matching, validity/noise or all-sector stability certificate. Author ordinary-argument and arithmetic/provenance review only.')
ans=dict(status='actual_fractional_comparison_and_exact_reference_residual_reviewed_actual_action_open',
 finished_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(raw),
 review_source_sha256=digest(Path(__file__)),review_source_utf8=Path(__file__).read_text(encoding='utf-8'),
 checked_inputs=len(d['input_sha256_before']),protected_paths_unchanged=len(baseline),
 cell_count=len(cells),metric_identity_fixtures=fixtures,review=reviews,
 constants={k:enc(v) for k,v in c.items()},
 exact_analytic_first_reference_screen=True,implemented_first_actual_inverse=False,
 complete_original_admission=False,goal_complete=False)
out.write_bytes((json.dumps(ans,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=ans['status'],checked_inputs=ans['checked_inputs'],cells=len(cells),
 reference_residual=c['complete_exact_analytic_first_reference_residual'].str(12),
 reference_energy_error=c['reference_inverse_energy_error'].str(12)),indent=2))
