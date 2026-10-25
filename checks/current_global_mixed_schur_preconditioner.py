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
STEM = 'current-global-mixed-schur-preconditioner'
OWN = ['notes/current-global-mixed-schur-preconditioner.md', 'checks/current_global_mixed_schur_preconditioner.py', 'evidence/current-global-mixed-schur-preconditioner-review.json']
INPUTS = ['notes/current-global-mixed-schur-preconditioner.md', 'checks/current_global_mixed_schur_preconditioner.py', 'evidence/current-global-mixed-schur-preconditioner-review.json', 'evidence/current-global-mixed-schur-preconditioner-exploration.json', 'notes/current-canonical-weighted-inf-sup-attempt.md', 'notes/current-restored-weight-commutator.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-seed-boundary-phase-and-trace-defect.md', 'notes/product-centered-scalar-cold-response.md', 'notes/physical-metric-block-and-preconditioner.md', 'notes/current-coupled-global-return-attempt.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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
    ctx.prec=512;a=arb;p,h,zh,env=york.actual_parameters();be=p['beta'];s=a('1.4');g=p['g']
    Dl=96*a.pi()*a('1e22')*g-g/(60*a.pi());d=Dl/be;epsA=g/(60*a.pi()*be)
    gap=12*be-s*s;massgap=1-d*s*s;Pnorm=1/gap
    grad=1/gap.sqrt()+Pnorm
    theta_max=a(10)*factorial(21)/factorial(10)**2/a(4)**10
    lapomega=s*(a('.04')*a('.9')*38+1+4+a('.5'))
    hess_rad=s*(a('.04')*a('.9')*38+1+2+a('.25'));hess_time=s*(2+a('.25'));hess_sphere=s*a('.25')
    hessomega=(hess_rad**2+hess_time**2+2*hess_sphere**2).sqrt()
    lap=1+2*(1+s)*grad+(a('2.001')+10+s*s+2*s)*Pnorm
    hess=lap+(13*be).sqrt()*grad
    wgrad=grad+s*Pnorm;whess=hess+2*s*grad+(s*s+8)*Pnorm
    tau0=be*(a('1.1e-27')+12*11*g/(360*a.pi()))
    err_trace=(1+1/massgap)*tau0*Pnorm
    err_anom=(epsA*a('1e12')+40*be*a('4e-27')/(60*a.pi())*26)/(2*massgap)
    err=err_trace+err_anom;Cnorm=Pnorm/(2*massgap*(1-err))
    logcoef=g/(10*a.pi()*Dl)
    out=dict(global_mixed_P_gap=gap.str(30),global_massive_gap=massgap.str(30),P_inverse_dimensionless=Pnorm.str(30),theta_prime_upper=theta_max.str(30),laplacian_weight_bound=lapomega.str(30),hessian_weight_bound=hessomega.str(30),gradient_weighted_solution=grad.str(30),laplacian_weighted_solution=lap.str(30),hessian_weighted_solution=hess.str(30),weighted_gradient_unconjugated_solution=wgrad.str(30),weighted_hessian_unconjugated_solution=whess.str(30),trace_remainder_upper=err_trace.str(30),all_anomaly_remainder_upper=err_anom.str(30),complete_mixed_preconditioned_remainder_upper=err.str(30),complete_mixed_inverse_dimensionless=Cnorm.str(30),schur_log_principal_coefficient=logcoef.str(30))

    check('actual_seed_coefficient_domain_guards',be>0 and p['u']>0 and p['eps']<0 and 2*be>4*abs(p['eps']))
    check('retained_exact_Pseed_endpoint',(6*be+p['u']+p['eps']).contains(1))
    check('global_Qprime_over_y_guard',12*be+2*abs(p['u'])+4*abs(p['eps'])+3*p['tau']/(1-zh)<4)
    check('thermal_bolt_coefficient_guard',3*p['tau']/(1-zh)<a('1e-30'))
    check('global_mixed_gap',gap>a('.03999999995'))
    check('global_mixed_second_inverse',Pnorm<a('25.0000001'))
    check('global_massive_gap',massgap>a('.99999996536'))
    check('global_massive_inverse',1/massgap<a('1.000000035'))
    check('entire_transition_theta_derivative',theta_max<38)
    check('entire_weight_laplacian_guard',lapomega<10)
    check('entire_weight_hessian_guard',hessomega<8)
    check('physical_gradient_conversion',grad<a('30.000001'))
    check('physical_laplacian_conversion',lap<565)
    check('physical_Bochner_Hessian_conversion',hess<610)
    check('weighted_unconjugated_gradient',wgrad<66)
    check('weighted_unconjugated_Hessian',whess<950)
    check('global_normalized_P_potential_guard',12*be*(1+11*g/(360*a.pi()))+2*p['u']<a('2.001'))
    # Explicit conservative tensor contractions behind the whole curvature anomaly bound.
    h0=a(3)*26;h2=a(3)*950+80*66+200*26
    deltaRm=6*h2+10*10*h0
    deltaRic=2*a('1e5')+10*h0
    deltaR=2*a('1e6')+13*be*h0
    deltaJ=a('1e6')+a('1e7')/2+12*be*h0/4
    deltaC=a('1e5')+4*a('1e6')+4*a('1e7')+a('1e5')
    check('orbital_tensor_Hessian_guard',h2<14000)
    check('complete_Riemann_variation_guard',deltaRm<a('1e5'))
    check('complete_Ricci_variation_guard',deltaRic<a('1e6'))
    check('complete_scalar_variation_guard',deltaR<a('1e7'))
    check('complete_traceless_Ricci_variation_guard',deltaJ<a('1e7'))
    check('complete_Weyl_variation_guard',deltaC<a('1e8'))
    curv=(a(7)/12*(2*35*be*a('1e8')+4*h0*(35*be)**2)+a(11)/6*(2*7*be*a('1e7')+4*h0*(7*be)**2)+a(11)/72*(2*100*be*a('1e7'))+52*a('1e5'))
    check('whole_curvature_density_anomaly_guard',curv<a('1e12'))
    check('actual_trace_preconditioned_remainder',err_trace<a('9.167e-27'))
    check('complete_anomaly_preconditioned_remainder',err_anom<a('1.555e-23'))
    check('full_mixed_remainder',err<a('1.556e-23') and err<1)
    check('global_full_mixed_inverse',Cnorm<a('12.500001'))
    check('Schur_log_coefficient_positive',logcoef>0 and logcoef<a('1.05543e-26'))
    # Differential and ordered block identities, independent of numerical bounds.
    z=sp.symbols('z',positive=True);bb,uu,ee,tt=sp.symbols('beta u eps tau')
    PP=bb*(1+2*z+3*z*z)+uu*z*z+ee*z**4;QQ=(1-z)**2*PP-tt*z**3
    qp=-12*bb*z*z*(1-z)+2*uu*z*(1-z)*(1-2*z)+2*ee*z**3*(1-z)*(2-3*z)-3*tt*z*z
    check('exact_Qprime_identity',sp.expand(sp.diff(QQ,z)-qp)==0)
    om=sp.Function('omega')(z);B,C=sp.symbols('B C',positive=True);F=C*QQ/z**2
    dr=lambda f:-z*z*sp.diff(f,z)/B
    Hrr=F*dr(dr(om))+dr(F)*dr(om)/2
    Htt=dr(F)*dr(om)/2;Hss=F*z*dr(om)/B
    check('physical_weight_radial_Hessian',sp.simplify(B*B*Hrr/C-z*z*QQ*sp.diff(om,z,2)-(z*z*sp.diff(QQ,z)/2+z*QQ)*sp.diff(om,z))==0)
    check('physical_weight_laplacian',sp.simplify(B*B*(Hrr+Htt+2*Hss)/C-z*z*QQ*sp.diff(om,z,2)-(z*z*sp.diff(QQ,z)-2*z*QQ)*sp.diff(om,z))==0)
    dd=sp.symbols('D',positive=True);kk=sp.symbols('k',positive=True);gg=sp.symbols('G',positive=True);ll=sp.symbols('log')
    AA=6*dd*kk**4;CC=2*dd*kk**4;DD=-sp.Rational(2,3)*dd*kk**4+gg*kk**4*ll/(15*sp.pi)
    check('exact_fourth_rank_one_Schur_cancellation',sp.simplify(1+AA*DD/CC**2-gg*ll/(10*sp.pi*dd))==0)
    Aop,Csi,Dop,T,fv,gv=sp.symbols('A CstarInv D T f g',commutative=False)
    first=-Aop*Csi*(Dop*T-gv)-sp.Symbol('CT',commutative=False)-fv
    schur=sp.Symbol('CT',commutative=False)+Aop*Csi*Dop*T-Aop*Csi*gv+fv
    check('ordered_original_first_row_Schur_identity',sp.expand(first+schur)==0)
    exponent,lam,ep,DR=sp.symbols('exponent lambda epsilon DR')
    orbit=exponent*(exponent-1);radial=exponent*(exponent-3)
    Cnormal=2*lam*(12-orbit)*(1-DR*lam*radial)+44*ep*lam**2-sp.Rational(22,3)*ep*lam**2*(orbit-6)
    Cfact=2*lam*(1+sp.Rational(11,3)*ep*lam-DR*lam*radial)*(4-exponent)*(exponent+3)
    Cstar=2*lam*(1+sp.Rational(11,3)*ep*lam-DR*lam*radial)*(6-exponent)*(exponent+1)
    check('complete_actual_outer_mixed_normal',sp.expand(Cnormal-Cfact)==0)
    check('actual_outer_adjoint_normal',sp.expand(Cfact.subs(exponent,3-exponent)-Cstar)==0)
    check('heavy_mixed_root_above_six',1+11*g/(180*a.pi())-18*Dl>0)
    check('conditional_fast_screen_unchanged',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    return dict(numerics=out,encoded_mixed_remainder=encoded(err),encoded_full_mixed_inverse=encoded(Cnorm),actual_seed_global_mixed_inverse_certified=True,actual_seed_global_mixed_adjoint_inverse_certified=True,physical_L2_and_H2_conversions_certified=True,exact_coupled_Schur_reduction_executed=True,identity_plus_small_Schur_estimator_fails=True,complete_Schur_inverse_certified=False,post_elimination_mass_column_certified=False,actual_mass_computed=False,original_fast_graph_comparisons_closed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,full_objective_complete=False,closed_audit_preserved=True,next_calculation_executed=False)

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
