"""Available contact budget and residual identities; no physical metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
import sympy as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_seed_radial_green_enclosure as prior
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-vacuum-contact-feasibility-closeout'
PREVIOUS='current-bolt-adjoint-domain'
OWN=['notes/'+STEM+'.md','checks/current_vacuum_contact_feasibility_closeout.py','evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; FIXTURES=[]; encoded=prior.encoded

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)
prior.check=check

def decode(x):
    z=x['decimal_enclosure']; e=int(z['exponent'])
    return arb(str(z['mid'])+'e'+str(e),str(z['rad'])+'e'+str(e))

def identities():
    b,A0,A2,B2,C,G,J,F,R2=sp.symbols('b A0 A2 B2 C G J F R2',nonzero=True)
    action=-(B2/C)*A0*sp.exp(2*b)/2
    row=sp.diff(action,b,2).subs(b,0)*(-2*G/B2)
    check('canonical_vacuum_Hessian_density',sp.simplify(row-4*G*A0/C)==0)
    check('angular_vacuum_coefficient',sp.simplify(row.subs(A0,J/(2*sp.pi*B2))-2*G*J/(sp.pi*C*B2))==0)
    fourth=(B2/C)/2*(-B2*F/(480*sp.pi))*(C/B2)**2
    check('independent_fourth_contact_volume_conversion',sp.simplify(fourth+C*F/(960*sp.pi))==0)
    em=130+192+sp.Rational(1176,12)+sp.Rational(3000,720)+sp.Rational(7056,30240)+sp.Rational(8*15984,720)
    check('retained_EM_arithmetic_inside_604',em<604 and 604<1000,em)
    # Independent exact finite-matrix identity, not the quantum Hessian.
    H=sp.Matrix([[-7,2],[2,5]]); r=sp.Matrix([3,-11]); ell=sp.Matrix([5,sp.Rational(5,2)])
    e=-H.inv()*r; psi=H.inv()*ell
    ea=e+sp.Matrix([sp.Rational(1,13),sp.Rational(-1,17)])
    pa=psi+sp.Matrix([sp.Rational(-2,19),sp.Rational(1,23)])
    d=H*ea+r; db=ell-H*pa
    approx=(ell.T*ea-pa.T*d)[0]; error=(ell.T*e)[0]-approx
    check('two_residual_corrected_observable_identity',sp.simplify(error+(db.T*H.inv()*d)[0])==0)
    dual=lambda v:v[0]**2/7+v[1]**2/5
    check('two_residual_unit_form_error_fixture',error**2<=dual(db)*dual(d))
    D=sp.symbols('D',positive=True); eta=sp.Rational(1667,10**12)
    check('half_budget_residual_product',sp.simplify(D*(sp.Rational(8335,10**13)/D)-eta/2)==0)
    FIXTURES.append(dict(fixture='exact_normalizations_and_two_residual_identity',actual_operator_test=False))

def budgets(old):
    p=prior.parameters(); q=arb(p['q']); k=arb(2000)
    C0=(1+12*arb.pi()*k**6).sqrt(); b02=k**10*(C0-1)/6
    L=(p['B2']/b02).log()/2
    check('unchanged_original_scale_guard',L.abs_upper()<arb('1e-8'),encoded(L))
    U=2*(604*q*q+(1+2*L.abs_upper())*(q*q/6-arb(1)/30))/(arb.pi()*p['C']*p['B2'])
    Uold=2*(1000*q*q+(1+2*L.abs_upper())*(q*q/6-arb(1)/30))/(arb.pi()*p['C']*p['B2'])
    value=arb('1.2e-8')*U; derivative=arb('7.7e-6')*U
    cold=arb('1e-6')*U; fraction=cold/arb('1.72e-15')
    check('vacuum_contact_value_budget',value<arb('1e-28'),encoded(value))
    check('vacuum_contact_derivative_budget',derivative<arb('1e-28'),encoded(derivative))
    check('vacuum_contact_cold_materiality',fraction<arb('1e-14'),encoded(fraction))
    check('older_1000_bound_already_suffices',Uold*arb('7.7e-6')<arb('1e-28'))
    NUMERICS.append(dict(kind='available_vacuum_contact_budget',log_mu0B=encoded(L),known_bound_coefficient_604=encoded(U),older_bound_coefficient_1000=encoded(Uold),interior_value_allowance=encoded(value),interior_derivative_allowance=encoded(derivative),cold_a_s_coefficient=encoded(cold),fraction_of_massive_allowance=encoded(fraction),angular_zeta_derivative_evaluated=False,complete_constant_row_evaluated=False,outer_contact_estimate=False))
    base=next(v for v in old['numerics'] if v['kind']=='actual_current_bolt_normalization')
    h=decode(base['h']); ell=(1+h)/h; eta=arb('1.667e-9'); tol=eta/ell
    check('decoded_bolt_interval_positive',h>0)
    check('physical_threshold_not_correction',tol>arb('8e-43') and tol<arb('9e-43'),encoded(tol))
    check('physical_tolerance_normalization_identity',(tol*ell/eta).contains(1))
    NUMERICS.append(dict(kind='physical_bolt_metric_combination_tolerance',h=encoded(h),physical_functional_factor=encoded(ell),target=encoded(eta),metric_combination_tolerance=encoded(tol),half_residual_product_budget=encoded(eta/2),actual_correction_computed=False))
    FIXTURES.append(dict(fixture='current_parameters_and_contact_materiality',full_solution_feasibility_test=False))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; deps=prior.dependency_identity(); error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_domain_execution_passed',old['status']=='passed' and old['physical_bolt_adjoint_source_identified'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        identities(); budgets(old)
    except Exception: error=traceback.format_exc()
    after={p:digest(p) for p in paths}; depafter=prior.dependency_identity()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),known_vacuum_contact_budget_substituted=status=='passed',two_residual_bolt_certificate_derived=status=='passed',angular_zeta_derivative_evaluated=False,full_bolt_adjoint_computed=False,full_defect_bolt_pairing_evaluated=False,full_adjoint_residual_norm_bounded=False,physical_bolt_functional_dual_norm_quantified=False,all_local_renormalized_moments_computed=False,coincident_subtraction_executed=False,differentiated_response_operator_bounded=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run); out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]: log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
    if status!='passed': raise SystemExit(1)

if __name__=='__main__': main()
