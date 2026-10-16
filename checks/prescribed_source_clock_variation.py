"""S134 fixed laboratory source versus fixed global profile; no Lean."""
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / 'evidence/prescribed-source-clock-variation.json'
LOG = ROOT / 'evidence/prescribed-source-clock-variation.log'
DPS = 60
INITIAL_STEPS = ['2e-4', '1e-4']
STEPS = ['2e-5', '1e-5']
TOLERANCES = dict(action_identity='1e-45', first_difference='1e-14',
                  second_difference='1e-12', nonzero_witness_minimum='1e-6',
                  retiming_difference='1e-16', inequality_rounding='1e-40')
INPUTS = ['notes/source-clock-allocation.md', 'notes/source-clock-window.md',
          'notes/source-clock-controller.md', 'notes/curved-fermion-boundary.md',
          'notes/boundary-vertex-matching.md', 'notes/cutoff-clock-contact-matching.md',
          'notes/source-clock-radiation.md', 'notes/prescribed-source-clock-variation.md',
          'checks/prescribed_source_clock_variation.py']


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def write_utf8_lf(path, value):
    with path.open('w', encoding='utf-8', newline='\n') as stream:
        stream.write(value)


def record(message):
    with LOG.open('a', encoding='utf-8', newline='\n') as stream:
        stream.write(message + '\n')
    print(message, flush=True)


def zero(expr):
    if s.simplify(expr) != 0:
        raise AssertionError(str(expr))


def symbolic_checks():
    u, parameter, other = s.symbols('u parameter other', real=True)
    A, omega = s.symbols('A omega', real=True)
    y = s.symbols('y', positive=True)
    r, rp, h0, hp = s.symbols('r rp h hp', real=True)
    # Coordinate representative is varied at fixed argument, then composed
    # with the moving clock. Both its value and its u-rate cancel exactly.
    delta_theta = -r*h0
    delta_r = -rp*h0-r*hp/y
    zero(delta_theta+r*h0)
    zero(hp*r+y*(delta_r+rp*h0))
    fixed_rate = s.symbols('fixed_rate', real=True)
    ys = y+parameter*hp
    fixed_density = A*ys**2*(fixed_rate/ys)**2
    zero(fixed_density-A*fixed_rate**2)
    zero(s.diff(fixed_density, parameter))
    zero(s.diff(fixed_density, parameter, 2))

    t = s.Function('t')(u)
    profile = s.Function('Theta')
    density = A*s.diff(t, u)**2*s.diff(profile(t), t)**2
    euler = s.diff(density, t)-s.diff(s.diff(density, s.diff(t, u)), u)
    expected = -2*A*(s.diff(profile(t), t)**2*s.diff(t, u, 2)
                     +s.diff(t, u)**2*s.diff(profile(t), t)*s.diff(profile(t), t, 2))
    zero(euler-expected)

    h, k = u**4*(1-u)**3, u**3*(1-u)**4
    for v in (h, k):
        for order in range(3):
            zero(s.diff(v, u, order).subs(u, 0))
            zero(s.diff(v, u, order).subs(u, 1))
    global_time = s.symbols('global_time', real=True)
    theta = omega*(10*global_time**3-15*global_time**4+6*global_time**5)
    rate = s.diff(theta, global_time)
    r_at_u = rate.subs(global_time, u)
    rp_at_u = s.diff(rate, global_time).subs(global_time, u)
    direct_first = 2*A*(r_at_u**2*s.diff(h, u)+r_at_u*rp_at_u*h)
    integrated_first = s.integrate(s.expand(direct_first), (u, 0, 1))
    integrated_euler = s.integrate(s.expand(-2*A*r_at_u*rp_at_u*h), (u, 0, 1))
    base = s.integrate(s.expand(A*r_at_u**2), (u, 0, 1))
    zero(integrated_first-10*A*omega**2/1001)
    zero(integrated_first-integrated_euler)
    zero(base-10*A*omega**2/7)
    # Independently differentiate the displaced density at the polynomial
    # witness; compare with the general first and mixed Hessian formulas.
    clock = u+parameter*h+other*k
    displaced = A*s.diff(clock, u)**2*rate.subs(global_time, clock)**2
    from_displacement = s.diff(displaced, parameter).subs({parameter:0, other:0})
    zero(from_displacement-direct_first)
    rpp_at_u = s.diff(rate, global_time, 2).subs(global_time, u)
    mixed = 2*A*(r_at_u**2*s.diff(h,u)*s.diff(k,u)
                  +2*r_at_u*rp_at_u*(h*s.diff(k,u)+k*s.diff(h,u))
                  +(rp_at_u**2+r_at_u*rpp_at_u)*h*k)
    from_mixed = s.diff(displaced, parameter, other).subs({parameter:0, other:0})
    zero(from_mixed-mixed)
    return dict(fixed_source_first_second='identically zero',
                fixed_global_profile_euler=str(expected),
                quintic_witness_action=str(base),
                quintic_witness_first_derivative=str(integrated_first),
                endpoint_variation_jets='h,k and their first two derivatives vanish at 0 and 1',
                monotonicity_neighborhood='|s|,|v|<=1; |hprime|,|kprime|<=1/4 implies y>=1/2')


def integrate(fn):
    return mp.quad(fn, [mp.mpf(0), mp.mpf('0.5'), mp.mpf(1)])


def direction(u):
    return u**4*(1-u)**3


def direction_prime(u):
    return 4*u**3*(1-u)**3-3*u**4*(1-u)**2


def other_direction(u):
    return u**3*(1-u)**4


def other_prime(u):
    return 3*u**2*(1-u)**4-4*u**3*(1-u)**3


def action_differences():
    A, omega = mp.mpf(7)/5, mp.mpf(2)/5
    first_tol = mp.mpf(TOLERANCES['first_difference'])
    second_tol = mp.mpf(TOLERANCES['second_difference'])
    rows = []
    # The second profile is not a returned switch; it supplies an independent
    # nonconstant-rate control on the same endpoint-preserving clock family.
    for name in ['quintic', 'quadratic']:
        if name == 'quintic':
            rate = lambda t: 30*omega*t*t*(1-t)**2
            rate_p = lambda t: 60*omega*t*(1-t)*(1-2*t)
            rate_pp = lambda t: 60*omega*(1-6*t+6*t*t)
        else:
            rate = lambda t: omega*t
            rate_p = lambda t: omega
            rate_pp = lambda t: mp.mpf(0)

        def action(a, b):
            def density(u):
                clock = u+a*direction(u)+b*other_direction(u)
                velocity = 1+a*direction_prime(u)+b*other_prime(u)
                if velocity <= 0:
                    raise AssertionError('nonmonotone displaced clock')
                return A*velocity**2*rate(clock)**2
            return integrate(density)

        first = integrate(lambda u: -2*A*rate(u)*rate_p(u)*direction(u))
        def hessian(h, hp, k, kp):
            return integrate(lambda u: 2*A*(rate(u)**2*hp(u)*kp(u)
                           +2*rate(u)*rate_p(u)*(h(u)*kp(u)+k(u)*hp(u))
                           +(rate_p(u)**2+rate(u)*rate_pp(u))*h(u)*k(u)))
        pure = hessian(direction, direction_prime, direction, direction_prime)
        mixed = hessian(direction, direction_prime, other_direction, other_prime)
        center = action(0, 0)
        if name == 'quintic':
            assert abs(center-10*A*omega**2/7) < mp.mpf(TOLERANCES['action_identity'])
            assert abs(first-10*A*omega**2/1001) < mp.mpf(TOLERANCES['action_identity'])
            assert abs(first) > mp.mpf(TOLERANCES['nonzero_witness_minimum'])
        # Retain the original resolution and verify the O(step^2) first-
        # difference convergence that motivated refining, rather than
        # relaxing the first-derivative tolerance after its initial failure.
        original_errors = []
        for raw_step in INITIAL_STEPS:
            step = mp.mpf(raw_step)
            error = abs((action(step,0)-action(-step,0))/(2*step)-first)
            original_errors.append(error)
        convergence_ratio = original_errors[0]/original_errors[1]
        assert abs(convergence_ratio-4) < mp.mpf('1e-6'), convergence_ratio
        differences = []
        for raw_step in STEPS:
            step = mp.mpf(raw_step)
            plus, minus = action(step, 0), action(-step, 0)
            fd1 = (plus-minus)/(2*step)
            fd2 = (plus-2*center+minus)/step**2
            fdm = (action(step,step)-action(step,-step)-action(-step,step)+action(-step,-step))/(4*step**2)
            errors = [abs(fd1-first), abs(fd2-pure), abs(fdm-mixed)]
            assert errors[0] < first_tol, errors
            assert max(errors[1:]) < second_tol, errors
            differences.append(dict(step=raw_step, first_error=str(errors[0]),
                                    pure_second_error=str(errors[1]), mixed_error=str(errors[2])))
        rows.append(dict(profile=name, A=str(A), omega=str(omega), action=str(center),
                         fixed_global_profile_first=str(first), pure_second=str(pure), mixed_second=str(mixed),
                         fixed_laboratory_source_first='0', fixed_laboratory_source_second='0',
                         original_coarse_steps=INITIAL_STEPS,
                         original_coarse_first_errors=list(map(str,original_errors)),
                         original_coarse_error_ratio=str(convergence_ratio),
                         differences=differences))
    return rows


def finite_pullback():
    # A closed-form finite clock and its inverse avoid solving or assuming
    # an inverse derivative. mp.diff differentiates the composed source
    # independently in t before the t-coordinate action is integrated.
    A, omega = mp.mpf(7)/5, mp.mpf(2)/5
    theta_b = lambda u: omega*(10*u**3-15*u**4+6*u**5)
    rate_b = lambda u: 30*omega*u**2*(1-u)**2
    base = integrate(lambda u: A*rate_b(u)**2)
    tolerance = mp.mpf(TOLERANCES['action_identity'])
    rows = []

    def transformed_actions(parameter):
        clock = lambda u: u/(1+parameter*(1-u))
        inverse = lambda t: (1+parameter)*t/(1+parameter*t)
        velocity = lambda u: (1+parameter)/(1+parameter*(1-u))**2
        global_profile = lambda t: theta_b(inverse(t))
        u_action = integrate(lambda u: A*velocity(u)**2*mp.diff(global_profile,clock(u))**2)
        t_action = integrate(lambda t: A*velocity(inverse(t))*mp.diff(global_profile,t)**2)
        return u_action, t_action, clock, inverse, global_profile

    for raw in ['-0.4','-0.15','0.2','0.5']:
        parameter = mp.mpf(raw)
        u_action, t_action, clock, inverse, profile = transformed_actions(parameter)
        source_error = max(abs(profile(clock(mp.mpf(i)/16))-theta_b(mp.mpf(i)/16)) for i in range(17))
        inverse_error = max(abs(inverse(clock(mp.mpf(i)/16))-mp.mpf(i)/16) for i in range(17))
        errors = [abs(u_action-base), abs(t_action-base), source_error, inverse_error]
        assert max(errors) < tolerance, errors
        rows.append(dict(parameter=raw, u_action=str(u_action), t_action=str(t_action),
                         u_action_error=str(errors[0]), t_action_error=str(errors[1]),
                         source_composition_error=str(source_error), inverse_error=str(inverse_error)))
    # Independently take finite differences of the compensated finite action.
    differences = []
    for raw_step in STEPS:
        step = mp.mpf(raw_step)
        plus = transformed_actions(step)[1]
        minus = transformed_actions(-step)[1]
        first, second = (plus-minus)/(2*step), (plus-2*base+minus)/step**2
        assert abs(first) < mp.mpf(TOLERANCES['first_difference'])
        assert abs(second) < mp.mpf(TOLERANCES['second_difference'])
        differences.append(dict(step=raw_step, first_error=str(abs(first)), second_error=str(abs(second))))
    return dict(base_action=str(base), finite_clocks=rows, compensated_action_differences=differences,
                coordinate_identity_does_not_make_physical_clock_embedding_pure_gauge=True)


def retiming():
    a = mp.mpf(1)/5
    sigma = mp.mpf(7)/10
    T = mp.mpf(1)
    tolerance = mp.mpf(TOLERANCES['retiming_difference'])
    allowance = mp.mpf(TOLERANCES['inequality_rounding'])
    rows = []
    for b in [mp.mpf(3)/10, -mp.mpf(1)/10, a]:
        # Both clocks are derived independently from their laboratory
        # velocities: t0(u)=(exp(a*u)-1)/a and t(u)=(exp(b*u)-1)/b.
        t0 = lambda u: mp.expm1(a*u)/a
        inverse = lambda t: mp.log1p(b*t)/b
        F = lambda t: mp.expm1(a*inverse(t))/a
        exact = lambda t: (1+b*t)**(a/b-1)
        theta_b = lambda u: mp.sin(sigma*t0(u))
        actual_source = lambda t: theta_b(inverse(t))
        m = min(mp.mpf(1), 1+b*T)
        L, delta, R0, R1 = abs(a), abs(b-a)*T, sigma, sigma**2
        M0 = max(1+a*T, 1+a*F(T))
        errors, estimates = [], []
        for i in range(21):
            t = T*i/20
            ode_rhs = (1+a*F(t))/(1+b*t)
            fd = (F(t+mp.mpf('1e-8'))-F(t-mp.mpf('1e-8')))/(2*mp.mpf('1e-8'))
            derivative_error = abs(fd-ode_rhs)
            assert derivative_error < tolerance, derivative_error
            assert abs(exact(t)-ode_rhs) < allowance
            source_exact = mp.sin(sigma*F(t))
            assert abs(actual_source(t)-source_exact) < allowance
            source_rate = mp.diff(actual_source,t)
            predicted_rate = sigma*mp.cos(sigma*F(t))*ode_rhs
            assert abs(source_rate-predicted_rate) < allowance
            E = delta*mp.expm1(L*t/m)/L
            clock_error = abs(F(t)-t)
            velocity_error = abs(ode_rhs-1)
            angle_error = abs(source_exact-mp.sin(sigma*t))
            rate_error = abs(source_rate-sigma*mp.cos(sigma*t))
            bounds = [E, (L*E+delta)/m, R0*E, R1*M0/m*E+R0*(L*E+delta)/m]
            for actual, bound in zip([clock_error,velocity_error,angle_error,rate_error],bounds):
                assert actual <= bound+allowance, (actual,bound)
            errors.append(derivative_error)
            estimates.append(dict(t=str(t), F=str(F(t)), clock_error=str(clock_error),
                                  angle_error=str(angle_error), rate_error=str(rate_error),
                                  clock_bound=str(E), rate_bound=str(bounds[-1])))
        rows.append(dict(reference_linear_slope=str(a), actual_linear_slope=str(b),
                         m=str(m), L=str(L), delta=str(delta), M0=str(M0), R0=str(R0), R1=str(R1),
                         maximum_ode_derivative_difference_error=str(max(errors)), estimates=estimates,
                         reference_trajectory_case=(a==b), affine_feedback_law_asserted=False))
    # Independent elementary constant-reference branch L=0 of Gronwall.
    for b in [mp.mpf(1)/4, -mp.mpf(1)/4]:
        m = min(mp.mpf(1),1+b*T)
        delta = abs(b)*T
        for i in range(21):
            t = T*i/20
            F = mp.log1p(b*t)/b
            assert abs(F-t) <= delta*t/m+allowance
    return dict(exact_open_loop_prescription='theta_b(u)=Theta0(t0(u)); F=t0(u(t)); Fprime=z0(F)/z(t)',
                cases=rows, constant_reference_L_zero_branch_checked=True,
                actual_controller_or_physical_state_implemented=False)


def main():
    mp.mp.dps = DPS
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    hashes = {name:sha(ROOT/name) for name in INPUTS}
    result = dict(story='S134', started_utc=started, dps=DPS, finite_difference_steps=STEPS,
                  initial_finite_difference_steps=INITIAL_STEPS,
                  refinement='steps divided by ten after initial first-difference truncation failure; no tolerance relaxed',
                  predeclared_tolerances=TOLERANCES, input_sha256=hashes,
                  hash_policy='SHA-256 after CRLF-to-LF normalization; no other byte changes',
                  execution=dict(executable=sys.executable, argv=sys.argv,
                                 interpreter_argv=sys.orig_argv, cwd=str(Path.cwd()),
                                 command=subprocess.list2cmdline(sys.orig_argv)),
                  runtime=dict(python=sys.version, platform=platform.platform(), sympy=s.__version__,mpmath=mp.__version__))
    if RECEIPT.exists():
        previous = json.loads(RECEIPT.read_text(encoding='utf-8'))
        if previous.get('status') == 'failed':
            record('Retained prior failed receipt: '+json.dumps(previous,ensure_ascii=False))
    record('S134 prescribed source clock variation run: '+started)
    record('Execution: '+json.dumps(result['execution'],ensure_ascii=False))
    record('Predeclared checks: '+json.dumps(dict(dps=DPS,steps=STEPS,tolerances=TOLERANCES)))
    try:
        result['symbolic'] = symbolic_checks()
        record('Ordinary identities and quintic counterexample: PASS')
        result['fixed_global_action_differences'] = action_differences()
        record('Independent displaced-action first/pure/mixed differences: PASS')
        result['finite_fixed_source_pullback'] = finite_pullback()
        record('Finite transformed-source pullback and compensated-action differences: PASS')
        result['laboratory_open_loop_retiming'] = retiming()
        record('Open-loop retiming ODE, source rate and Gronwall inequalities: PASS')
        assert hashes == {name:sha(ROOT/name) for name in INPUTS}, 'input changed during execution'
        result['status'] = 'pass_fixed_laboratory_contact_zero_clock_force_changed_profile_protocol'
        result['claims'] = dict(fixed_source_intrinsic_direct_clock_force_zero=True,
                               varying_clock_alone_selects_rule_D=False,
                               rule_D_oscillator_constraint_derived=False,
                               zeta_s_assigned_or_set_to_zero=False,
                               extrinsic_contact_removed=False,
                               quantum_or_field_dependent_terms_removed=False,
                               complete_controller_constraints_derived=False,
                               physical_source_implementation_proved=False,
                               s074_gate_closed=False, area_budget_evaluated=False,
                               physical_opening_achieved=False, lean_invoked=False)
        first_errors = [mp.mpf(d['first_error']) for row in result['fixed_global_action_differences'] for d in row['differences']]
        second_errors = [mp.mpf(d[key]) for row in result['fixed_global_action_differences'] for d in row['differences'] for key in ['pure_second_error','mixed_error']]
        result['summary'] = dict(action_profiles=2, finite_pullback_clocks=4, retiming_cases=3,
                               maximum_first_difference_error=str(max(first_errors)),
                               maximum_second_difference_error=str(max(second_errors)))
        write_utf8_lf(RECEIPT,json.dumps(result,indent=2,ensure_ascii=False)+'\n')
        record('Summary: '+json.dumps(result['summary']))
        record('Normalized input hashes: '+json.dumps(hashes))
        record('Terminal status: PASS, exit 0')
    except Exception:
        result['status'] = 'failed'
        result['failure_traceback'] = traceback.format_exc()
        write_utf8_lf(RECEIPT,json.dumps(result,indent=2,ensure_ascii=False)+'\n')
        record(result['failure_traceback'].rstrip())
        record('Failure input hashes: '+json.dumps(hashes))
        record('Terminal status: FAIL, exit 1; failure retained in this append-only log')
        raise


if __name__ == '__main__':
    main()
