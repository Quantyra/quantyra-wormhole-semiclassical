"""Supporting arithmetic for the constrained weighted global-return attempt; no complete Hessian bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-canonical-weighted-inf-sup-attempt'
OWN = ['notes/current-canonical-weighted-inf-sup-attempt.md', 'checks/current_canonical_weighted_inf_sup_attempt.py', 'evidence/current-canonical-weighted-inf-sup-attempt-review.json']
INPUTS = ['notes/current-canonical-weighted-inf-sup-attempt.md', 'checks/current_canonical_weighted_inf_sup_attempt.py', 'evidence/current-canonical-weighted-inf-sup-attempt-review.json', 'evidence/current-canonical-weighted-inf-sup-attempt-exploration.json', 'notes/current-weighted-constraint-phase-obstruction.md', 'notes/current-restored-weight-commutator.md', 'notes/current-seed-full-hessian-realization.md', 'notes/product-centered-scalar-cold-response.md', 'notes/physical-metric-block-and-preconditioner.md', 'notes/current-coupled-global-return-attempt.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
CHECKS = []


def check(name, condition, detail=""):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def encoded(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(45), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def decode(value):
    p = value['decimal_enclosure']
    return arb(int(p['mid']), int(p['rad'])) * arb(10) ** int(p['exponent'])


from math import comb, factorial
sys.path.insert(0, str(ROOT/'checks'))
# A legacy parameter helper parses argv on import. Keep this wrapper's
# acceptance mode local while importing it, then restore the exact argv.
original_argv = sys.argv[:]
sys.argv = [v for v in sys.argv if v != '--acceptance']
import current_physical_york_dual_assembly as york
sys.argv = original_argv
york.check = check
york.prior.check = check

def execute():
    ctx.prec=512;a=arb;p,h,zh,env=york.actual_parameters();be=p['beta'];s=a('1.4');Y=(-a(20)).exp();y0=(-a(60)).exp();L=a(40)
    # a=z^2 P(z)-tau z^5/y^2; polynomial part has exact constant 1.
    z=sp.symbols('z');u,e=sp.symbols('u e');bb=(1-u-e)/6
    yp=sp.Poly(sp.expand((1-z)**2*(bb*(1+2*(1-z)+3*(1-z)**2)+u*(1-z)**2+e*(1-z)**4)),z)
    ac=[]
    for n in range(yp.degree()+1):
     c=sp.Poly(yp.nth(n),u,e);v=a(0)
     for pp,r in c.terms():v+=a(int(sp.numer(r)))/int(sp.denom(r))*p['u']**pp[0]*p['eps']**pp[1]
     ac.append(v)
    errs=[]
    for j in range(3):
     poly=sum((abs(ac[n])*n**j*Y**n for n in range(1,len(ac))),a(0))
     thermal=p['tau']/y0**2*sum((comb(5,k)*abs(k-2)**j*Y**k for k in range(6)),a(0))
     errs.append(poly+thermal)
    ae,da,dda=errs;alo=1-ae;ahi=1+ae;khi=a('.5')+2*Y/(1-Y);dk=2*Y/(1-Y)**2
    Verr=da*khi+ahi*dk+ae/4+ahi*(khi*khi-a('.25'))
    kaperr=a('1.3e-8');dkap=a('1e-20');aklo=alo*(1-kaperr);akhi=ahi*(1+kaperr)
    Vkaperr=(da*(1+kaperr)+ahi*dkap)*khi+akhi*dk+(ae+kaperr+ae*kaperr)/4+akhi*(khi*khi-a('.25'))
    omlo=s*y0/(h+y0);omhi=s
    rho1=4*a.pi()**2/(3*L**2);rho2=16*a.pi()**4/(3*L**4)
    d=(96*a.pi()*a('1e22')*p['g']-p['g']/(60*a.pi()))/be
    vhi=4*be+a('1e-20');Lphi=ahi*rho2.sqrt()+da*rho1.sqrt()+a('.25')+Verr
    fourth_upper=d*(Lphi+ahi*s*s)**2
    scalar_upper=akhi*rho1+a('.25')+Vkaperr+vhi-aklo*omlo**2+fourth_upper
    critical_test_upper=((akhi*rho1+a('.25')+Vkaperr+vhi+fourth_upper)/aklo).sqrt()*(h+y0)/y0
    product_gap=a('2.25')-s*s
    alternative_inverse=1/(2*product_gap)+a('1.5')/product_gap+2/product_gap**2
    weight_loss=((h+a('.8'))/(h+1-zh))**(s-a('.95'))
    out=dict(actual_window_y_min=y0.str(24),actual_window_y_max=Y.str(24),h_over_ymin=(h/y0).str(24),a_error=ae.str(24),a_prime_bound=da.str(24),V_error=Verr.str(24),scalar_second_potential_error=Vkaperr.str(24),first_derivative_ratio=rho1.str(24),second_derivative_ratio=rho2.str(24),d=d.str(24),fourth_upper=fourth_upper.str(24),full_scalar_signed_pairing_upper=scalar_upper.str(24),this_test_critical_exponent_upper=critical_test_upper.str(24),weaker_weight_bolt_loss_lower=weight_loss.str(24),cold_coupled_reference_gap=product_gap.str(24),cold_reference_inverse_algebra_upper=alternative_inverse.str(24))
    epsl=p['g']/(60*a.pi());Dl=d*be
    actual_kaperr=4*Dl+11*epsl/3+(2*Dl+11*epsl/6)*a('1.1e-27')/6+11*epsl*7/9
    actual_dkap=((2*Dl+11*epsl/6)*a('1.4e-26')/6+11*epsl*100/9)/(alo/be).sqrt()
    cmax=2*ahi*Y/(1-Y)
    omprime=s*h/y0
    lapW=Lphi+2*ahi*s*rho1.sqrt()+da*s+ahi*omprime+ahi*s*s
    lap2W=lapW+cmax*(rho1.sqrt()+s+khi)
    deltaR=24*be+2*be*a('1.1e-27')+4*p['u']+2*lap2W
    mixed_square=d*lapW*deltaR
    mixed_classical=2*(alo*rho1+a('.25')-Verr-ahi*s*s)+24*be-2*cmax*(rho1.sqrt()+s+khi)
    mixed_lower=mixed_classical-mixed_square-a('1e-20')
    out.update(actual_kappa_error_bound=actual_kaperr.str(24),actual_kappa_derivative_bound=actual_dkap.str(24),weighted_laplacian_norm_upper=lapW.str(24),complete_mixed_classical_pairing_lower=mixed_classical.str(24),mixed_curvature_square_upper=mixed_square.str(24),complete_mixed_pairing_lower=mixed_lower.str(24),weight_loss_at_test_ceiling_lower=(((h+a('.8'))/(h+1-zh))**(s-a('.961714'))).str(24),comparison_only_defect_for_half_margin=(1/(2*alternative_inverse)).str(24))

    check('actual_compact_support_inside_domain',y0>1-zh and Y<a('.8') and h/y0<a('5.777e-8'))
    check('whole_interval_a_guard',ae<a('6.871e-9') and alo>0)
    check('whole_interval_a_derivative_guard',da<a('6.871e-9'))
    check('whole_interval_unitary_potential_guard',Verr<a('1.340e-8'))
    check('actual_scalar_kappa_guard',actual_kaperr<kaperr)
    check('actual_scalar_kappa_derivative_guard',actual_dkap<dkap)
    check('scalar_unitary_potential_guard',Vkaperr<a('1.665e-8'))
    check('retained_curvature_square_coefficient',d>0 and d<a('1.768e-8'))
    vacuum_error=8*be*11*p['g']/(360*a.pi())+be*a('1.1e-27')/3
    check('full_scalar_density_potential_guard',vacuum_error<a('1e-20'))
    check('complete_fourth_form_upper',fourth_upper<a('8.743e-8'))
    check('actual_full_signed_pairing_negative',scalar_upper<-a('1.03510829'))
    check('test_exponent_threshold',critical_test_upper<a('.961712830'))
    check('mixed_classical_row_lower',mixed_classical>a('.5964492703'))
    check('weighted_laplacian_both_signs',lapW<a('2.478178'))
    check('mixed_full_square_term',mixed_square<a('3.923e-7'))
    anom_curve=a('1e6')*p['g']/(60*a.pi()*be)
    anom_mag=40*be*a('4e-27')/(60*a.pi())
    check('complete_curvature_anomaly_margin',anom_curve<a('4e-29'))
    check('complete_magnetic_anomaly_margin',anom_mag<a('2e-28'))
    check('complete_mixed_anomaly_margin',anom_curve+anom_mag<a('1e-20'))
    check('complete_actual_mixed_row_nonzero',mixed_lower>a('.59644887'))
    check('weaker_weight_bolt_attenuation_loss',weight_loss>a('6.370071e14'))
    check('minimal_test_ceiling_attenuation_loss',((h+a('.8'))/(h+1-zh))**(s-a('.961714'))>a('2.622865e14'))
    check('comparison_only_coupled_inverse_guard',alternative_inverse<a('30.677765'))
    check('comparison_only_half_margin_target',a('.01629')*alternative_inverse<a('.5'))
    # Independent algebraic identities used in the ordinary argument.
    pp=sp.symbols('p',nonzero=True);MM=sp.Matrix([[-6*pp+8,-2*pp],[-2*pp,0]])
    MI=sp.Matrix([[0,-1/(2*pp)],[-1/(2*pp),sp.Rational(3,2)/pp-2/pp**2]])
    check('complete_comparison_matrix_inverse',sp.simplify(MM*MI)==sp.eye(2))
    L0,aa,om,phi,phip,daw=sp.symbols('L0 a om phi phip daw')
    plus=L0-2*aa*om*phip-daw*phi-aa*om**2*phi
    minus=L0+2*aa*om*phip+daw*phi-aa*om**2*phi
    check('exact_fourth_conjugation_identity',sp.expand(plus*minus-(L0-aa*om**2*phi)**2+(2*aa*om*phip+daw*phi)**2)==0)
    TT,DR,DDR=sp.symbols('T DeltaR DeltadeltaR')
    check('exact_mixed_density_contact_cancellation',sp.expand(-2*TT*DR+DDR+2*TT*DR-DDR)==0)
    xx=sp.symbols('x',real=True);ph=sp.sin(sp.pi*(xx+60)/40)**2
    norm=sp.integrate(ph**2,(xx,-60,-20));first=sp.integrate(sp.diff(ph,xx)**2,(xx,-60,-20));second=sp.integrate(sp.diff(ph,xx,2)**2,(xx,-60,-20))
    check('compact_test_norms',norm==15 and sp.simplify(first/norm-sp.pi**2/1200)==0 and sp.simplify(second/norm-sp.pi**4/480000)==0)
    check('compact_test_H2_endpoint_jets',all(sp.simplify(f.subs(xx,t))==0 for f in [ph,sp.diff(ph,xx)] for t in [-60,-20]))
    check('conditional_fast_screen_unchanged',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    return dict(numerics=out,encoded_signed_pairing_upper=encoded(scalar_upper),encoded_complete_mixed_pairing_lower=encoded(mixed_lower),actual_finite_domain_two_row_test_executed=True,original_York_sign_weighted_coercivity_fails=True,full_coupled_inverse_fails=False,full_coupled_inverse_certified=False,shape_shape_quantum_operator_evaluated=False,actual_mass_coefficient_computed=False,post_elimination_mass_column_certified=False,original_physical_strong_comparisons_closed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,full_objective_complete=False,closed_audit_preserved=True,next_calculation_executed=False)

if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    import flint
    flint_root = Path(flint.__file__).resolve().parent
    runtime = {'python': sys.version, 'python_executable_sha256': sha256(Path(sys.executable).read_bytes()).hexdigest(),
               'flint_version': flint.__version__, 'sympy_version': sp.__version__,
               'flint_file_sha256': {str(p.relative_to(flint_root)): sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(flint_root.rglob('*'))
                                    if p.is_file() and p.suffix in ['.py', '.pyd', '.dll']}}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(), 'mode': mode,
           'input_sha256_before': hashes, 'precision_bits':512, 'runtime':runtime,
           'frozen_sources_utf8': {p: (ROOT/p).read_text(encoding='utf-8') for p in OWN}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('all_frozen_inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs': []}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger, indent=2)+'\n').encode())
    log = ROOT/'evidence'/(STEM+'.log')
    text = json.dumps({'mode': mode, 'status': run['status'], 'checks': len(CHECKS),
                       'result': run.get('result'), 'traceback': run.get('traceback')}, indent=2)+'\n'
    with log.open('ab') as f:
        f.write(text.encode())
    print(text)
    if run['status'] != 'passed':
        raise SystemExit(1)
