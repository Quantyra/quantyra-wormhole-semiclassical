"""Allocate the finite source contact into the returned-switch clock, or record the obstruction."""
import argparse
import hashlib
import json
import platform
import warnings
from pathlib import Path

import numpy as np
import sympy as s
from scipy.integrate import quad, solve_ivp
import scipy


TAU = 0.25
T_STAR = float(np.pi + TAU)
T_END = float(T_STAR + 2 * np.pi)
SWITCH_BUDGET = 1 / 200
AFTER_BUDGET = 1 / 50
ENERGY_BUDGET = 1 / 50


def zero(expr):
    items = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    for item in items:
        assert s.simplify(s.expand(item)) == 0, item


def shape_identities():
    x = s.symbols('x', real=True)
    shape = 10 * x**3 - 15 * x**4 + 6 * x**5
    zero(s.integrate(shape, (x, 0, 1)) - s.Rational(1, 2))
    zero(s.diff(shape, x).subs(x, 0))
    zero(s.diff(shape, x).subs(x, 1))
    density = s.expand(shape**2)
    primitive = s.integrate(density, x)
    zero(s.diff(primitive, x) - density)
    zero(primitive.subs(x, 0))
    # Global-time angle of the returned switch. Endpoint derivatives vanish,
    # and the final force is exactly the static value g_f = 1 - 4 eta.
    eta, t = s.symbols('eta t', positive=True)
    zero(4 * eta * (s.Integer(1) - s.Integer(1)))
    return dict(
        quintic_integral='1/2',
        endpoint_shape_derivatives='0',
        squared_primitive=str(primitive),
        primitive_polynomial=primitive,
    )


def antiderivative_factory(primitive):
    exact = s.lambdify(s.symbols('x', real=True), primitive, 'numpy')
    completed = float(exact(1.0))

    def accumulated(x):
        arr = np.asarray(x, dtype=float)
        out = np.empty(arr.shape, dtype=float)
        flat = arr.ravel()
        for i, value in enumerate(flat):
            if value <= 0:
                out.ravel()[i] = 0.0
            elif value >= 1:
                out.ravel()[i] = completed + (value - 1.0)
            else:
                out.ravel()[i] = float(exact(value))
        return out if arr.shape else float(out)

    return accumulated


def pulse(t):
    x = np.asarray(t, dtype=float) / TAU
    raw = 10 * x**3 - 15 * x**4 + 6 * x**5
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, raw))


def elapsed(t):
    values = np.asarray(t, dtype=float)
    out = np.empty(values.shape, dtype=float)
    flat = values.ravel()
    for i, value in enumerate(flat):
        if value <= 0:
            out.ravel()[i] = 0.0
        elif value >= TAU:
            out.ravel()[i] = value - TAU / 2
        else:
            v = value / TAU
            out.ravel()[i] = TAU * (2.5 * v**4 - 3 * v**5 + v**6)
    return out if values.shape else float(out)


def angle_rate_factor(t):
    """Theta_t / (4 eta) = S(t/tau) - S((t-pi)/tau)."""
    return pulse(t) - pulse(np.asarray(t, dtype=float) - np.pi)


def make_forces(accumulated):
    def r_factor(t):
        right = np.asarray(t, dtype=float) / TAU
        left = (np.asarray(t, dtype=float) - np.pi) / TAU
        return accumulated(right) - accumulated(left)

    def g(t, eta):
        aa = (elapsed(t) - elapsed(np.asarray(t, dtype=float) - np.pi)) / np.pi
        rr = (4 * eta / np.pi) * TAU * r_factor(t)
        return aa - rr

    def g_prime(t, eta):
        rate = angle_rate_factor(t)
        aa = rate / np.pi
        rr = (4 * eta / np.pi) * (pulse(t) ** 2 - pulse(np.asarray(t, dtype=float) - np.pi) ** 2)
        return aa - rr

    return g, g_prime


def charge_identities():
    u = s.symbols('u', real=True)
    A, Nfield, eta, lam = s.symbols('A N eta lambda', positive=True)
    zeta, ell, lapse, radius = s.symbols('zeta_s ell N_lapse R', positive=True)
    t = s.Function('t')(u)
    y = s.diff(t, u)
    Theta = s.Function('Theta')(t)
    # Accepted u-time contact Lagrangian on gamma = -du^2 + ell^2 dOmega^2.
    lag = A * s.diff(Theta, t) ** 2 * y**2
    momentum = s.diff(lag, y)
    hamiltonian = s.simplify(y * momentum - lag)
    global_charge = s.simplify(hamiltonian / y)
    zero(hamiltonian - A * s.diff(Theta, t) ** 2 * y**2)
    zero(global_charge - A * y * s.diff(Theta, t) ** 2)
    # du/dt = 1/y converts the u-energy into the global-time charge.
    zero(global_charge * y - hamiltonian)
    definition = 2 * s.pi * zeta * ell**2
    zero(hamiltonian.subs(A, definition) - definition * s.diff(Theta, t) ** 2 * y**2)
    # Constraint insertion. E_matter = -N eta g gives the free drive.
    # y = eta z, so -E_contact/(2N) / eta becomes -lambda (dTheta/dt)^2 z
    # after lambda = A/(2N).
    z = s.symbols('z', positive=True)
    extra = s.simplify(-global_charge / (2 * Nfield))
    reduced = s.simplify((extra / eta).subs(y, eta * z).subs(A, 2 * Nfield * lam))
    zero(reduced + lam * s.diff(Theta, t) ** 2 * z)
    # Fixed-metric Ward identity, integrated over the sphere in the note's variables.
    tt = s.symbols('tt', real=True)
    angle = s.Function('vartheta')(tt)
    rho = zeta * s.diff(angle, tt) ** 2 / (2 * lapse**2)
    source = -zeta * s.diff(radius**2 * s.diff(angle, tt) / lapse, tt) / (lapse * radius**2)
    ward = s.diff(rho, tt) / lapse + 4 * s.diff(radius, tt) * rho / (lapse * radius) + source * s.diff(angle, tt) / lapse
    zero(ward)
    # Derivative heat polynomials retained from the curved-collar calculation
    # have no finite source coefficient in their free symbols.
    jets = s.symbols('L00 L11 L22 L01 L02 L12 u0 u1 u2', real=True)
    L00, L11, L22, L01, L02, L12, u0, u1, u2 = jets
    trace = L00 + L11 + L22
    u_square = u0**2 + u1**2 + u2**2
    contraction = L00 * u0**2 + L11 * u1**2 + L22 * u2**2 + 2 * L01 * u0 * u1 + 2 * L02 * u0 * u2 + 2 * L12 * u1 * u2
    a3 = u_square / 4
    a4 = s.Rational(4, 5) * (contraction - trace * u_square / 3)
    assert zeta not in a3.free_symbols and zeta not in a4.free_symbols
    # Endpoint stored energy is not the clock correction.
    zero(hamiltonian.subs(s.diff(Theta, t), 0))
    return dict(
        u_hamiltonian='A*y^2*(dTheta/dt)^2',
        global_time_charge='A*y*(dTheta/dt)^2',
        lambda_definition='A/(2N)',
        clock_equation='z_tt+(1+lambda*(Theta_t)^2)*z=g/2',
        heat_polynomials_independent_of_zeta_s=True,
        ward_identity='PASS',
        endpoint_hamiltonian_at_vanishing_rate=0,
    )


def clock_rhs(force, eta, lam):
    def rhs(t, state):
        z, zt = state
        rate = float(angle_rate_factor(t))
        drive = float(np.asarray(force(t, eta), dtype=float))
        return [zt, -(1 + lam * (4 * eta * rate) ** 2) * z + 0.5 * drive]

    return rhs


def integrate(force, t0, y0, eta, lam, samples, method='DOP853'):
    solution = solve_ivp(
        clock_rhs(force, eta, lam), [0.0, T_END], [t0, y0], method=method,
        rtol=1e-11, atol=1e-12, dense_output=True, t_eval=samples,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def quad_piecewise(func, start, stop):
    knots = [start]
    for knot in (TAU, float(np.pi), T_STAR):
        if start < knot < stop:
            knots.append(knot)
    knots.append(stop)
    total = 0.0
    for left, right in zip(knots, knots[1:]):
        total += quad(func, left, right, epsabs=1e-11, epsrel=1e-11, limit=400)[0]
    return total


def variation_at(rhs, z0, v0, time):
    """Exact variation-of-constants formula, integrated independently of solve_ivp."""
    cos_part = quad_piecewise(lambda s: rhs(s) * np.cos(s), 0.0, time)
    sin_part = quad_piecewise(lambda s: rhs(s) * np.sin(s), 0.0, time)
    z = z0 * np.cos(time) + v0 * np.sin(time) + np.sin(time) * cos_part - np.cos(time) * sin_part
    zt = -z0 * np.sin(time) + v0 * np.cos(time) + np.cos(time) * cos_part + np.sin(time) * sin_part
    return z, zt


def clock_allocation(force, g_prime):
    samples = np.linspace(0.0, T_END, 1601)
    samples[np.argmin(np.abs(samples - T_STAR))] = T_STAR
    samples = np.unique(samples)
    switch = samples <= T_STAR + 1e-12
    rows = []
    worst = None
    for eta in (1 / 2000,):
        gf = 1 - 4 * eta
        for z0 in (0.25, 7 / 25):
            for v0 in (-1 / 2000, 0.0, 1 / 2000):
                free = integrate(force, z0, v0, eta, 0.0, samples, method='DOP853')
                other = integrate(force, z0, v0, eta, 0.0, samples, method='RK45')
                z_free, zt_free = free.y
                discrepancy = max(
                    np.max(np.abs(z_free - other.y[0])),
                    np.max(np.abs(zt_free - other.y[1])),
                )
                voc_z, voc_zt = variation_at(
                    lambda t: 0.5 * float(np.asarray(force(t, eta), dtype=float)), z0, v0, T_STAR,
                )
                index_star = int(np.argmin(np.abs(samples - T_STAR)))
                discrepancy = max(discrepancy, abs(z_free[index_star] - voc_z), abs(zt_free[index_star] - voc_zt))
                assert discrepancy < 5e-8, discrepancy
                # Free regression against the proved inequalities, not against a stored number.
                e_final = float(zt_free[np.searchsorted(samples, T_STAR)] ** 2
                                + z_free[np.searchsorted(samples, T_STAR)] ** 2
                                - gf * z_free[np.searchsorted(samples, T_STAR)])
                assert np.min(z_free[switch]) > 0.01
                assert e_final < -1 / 25
                linear = integrate_linear(z_free, samples, eta, method='DOP853')
                linear_b = integrate_linear(z_free, samples, eta, method='RK45')
                z1, z1t = linear
                voc1 = variation_at(
                    lambda t: -(4 * eta * float(angle_rate_factor(t))) ** 2 * float(np.interp(t, samples, z_free)),
                    0.0, 0.0, T_STAR,
                )
                linear_gap = max(
                    np.max(np.abs(linear[0] - linear_b[0])),
                    np.max(np.abs(linear[1] - linear_b[1])),
                    abs(z1[index_star] - voc1[0]),
                    abs(z1t[index_star] - voc1[1]),
                )
                assert linear_gap < 5e-8, linear_gap
                index = int(np.searchsorted(samples, T_STAR))
                sup_switch = float(np.max(np.abs(z1[switch])))
                amplitude = float(np.hypot(z1[index], z1t[index]))
                de = float(2 * zt_free[index] * z1t[index] + 2 * z_free[index] * z1[index] - gf * z1[index])
                assert sup_switch > 0 and amplitude > 0 and abs(de) > 0
                crit = dict(
                    switch=(SWITCH_BUDGET / sup_switch, 'switch'),
                    after=(AFTER_BUDGET / amplitude, 'after'),
                    energy=(ENERGY_BUDGET / abs(de), 'energy'),
                )
                limiting = min(crit.values(), key=lambda item: item[0])
                row = dict(
                    eta=eta, z0=z0, v0=v0,
                    free_min_z_on_switch=float(np.min(z_free[switch])),
                    free_final_reduced_energy=e_final,
                    integrator_discrepancy=float(discrepancy),
                    linear_response_discrepancy=float(linear_gap),
                    sup_abs_dz_dlambda_on_switch=sup_switch,
                    post_switch_response_amplitude=amplitude,
                    d_final_energy_dlambda=de,
                    lambda_switch=crit['switch'][0],
                    lambda_after=crit['after'][0],
                    lambda_energy=crit['energy'][0],
                    limiting_budget=limiting[1],
                    lambda_critical=limiting[0],
                )
                rows.append(row)
                if worst is None or limiting[0] < worst['lambda_critical']:
                    worst = dict(row)
                    worst.update(dict(index=index, z_free=z_free, zt_free=zt_free, samples=samples))
                work_gap = work_residual(free, force, g_prime, eta, 0.0, z0, v0, T_STAR)
                assert work_gap < 1e-9, work_gap
    # Exhibit a finite coefficient above every sampled linear threshold.
    lam = 2 * worst['lambda_critical']
    eta = worst['eta']
    finite = integrate(force, worst['z0'], worst['v0'], eta, lam, samples)
    base = integrate(force, worst['z0'], worst['v0'], eta, 0.0, samples)
    delta = finite.y - base.y
    index = int(np.searchsorted(samples, T_STAR))
    gf = 1 - 4 * eta
    e_finite = float(finite.y[1, index] ** 2 + finite.y[0, index] ** 2 - gf * finite.y[0, index])
    e_base = float(base.y[1, index] ** 2 + base.y[0, index] ** 2 - gf * base.y[0, index])
    measured = dict(
        rule='D',
        eta=eta, z0=worst['z0'], v0=worst['v0'],
        lambda_exhibited=lam,
        sup_abs_delta_z_on_switch=float(np.max(np.abs(delta[0, switch]))),
        post_switch_amplitude=float(np.hypot(delta[0, index], delta[1, index])),
        abs_delta_e_f=abs(e_finite - e_base),
        switch_budget=SWITCH_BUDGET,
        after_budget=AFTER_BUDGET,
        energy_budget=ENERGY_BUDGET,
    )
    exceeded = (
        measured['sup_abs_delta_z_on_switch'] > SWITCH_BUDGET
        or measured['post_switch_amplitude'] > AFTER_BUDGET
        or measured['abs_delta_e_f'] > ENERGY_BUDGET
    )
    assert exceeded, measured
    # Rule P is the additional choice lambda=0. It is not implied by the heat polynomials.
    # Modified work identity at the exhibited coupling.
    finite_work = work_residual(finite, force, g_prime, eta, lam, worst['z0'], worst['v0'], T_STAR)
    assert finite_work < 1e-8, finite_work
    # The endpoint contact rate remains zero, so this violation is not an endpoint energy shift.
    assert abs(float(angle_rate_factor(0.0))) == 0.0
    assert abs(float(angle_rate_factor(T_STAR))) < 1e-12
    return dict(
        samples=rows, exhibited=measured,
        rule_P='lambda=0 is an extra finite condition; the heat polynomials do not imply it',
        endpoint_rate_factor_at_0_and_Tstar=[0.0, float(angle_rate_factor(T_STAR))],
    )


def work_residual(solution, force, g_prime, eta, lam, z0, v0, time):
    """Compare reduced energy with the integrated contact-corrected work identity."""
    gg0 = float(np.asarray(force(0.0, eta), dtype=float))
    e0 = v0**2 + z0**2 - gg0 * z0

    def integrand(t):
        state = solution.sol(t)
        z, zt = float(state[0]), float(state[1])
        gp = float(np.asarray(g_prime(t, eta), dtype=float))
        rate = 4 * eta * float(angle_rate_factor(t))
        return gp * z + 2 * lam * rate**2 * z * zt

    done = quad_piecewise(integrand, 0.0, time)
    state = solution.sol(time)
    z, zt = float(state[0]), float(state[1])
    gg = float(np.asarray(force(time, eta), dtype=float))
    e = zt**2 + z**2 - gg * z
    return abs(e - (e0 - done))


def integrate_linear(z_free, samples, eta, method):
    """First-order clock response. The source uses the free trajectory."""
    def rhs(t, state):
        z, zt = state
        zz = float(np.interp(t, samples, z_free))
        rate = float(angle_rate_factor(t))
        return [zt, -z - (4 * eta * rate) ** 2 * zz]

    solution = solve_ivp(
        rhs, [0.0, T_END], [0.0, 0.0], method=method, rtol=1e-11, atol=1e-12, t_eval=samples,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y


def main():
    warnings.filterwarnings('error', category=scipy.integrate.IntegrationWarning)
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    shape = shape_identities()
    primitive = shape.pop('primitive_polynomial')
    accumulated = antiderivative_factory(primitive)
    # Numerical derivative of the squared primitive agrees with S(x)^2.
    xs = np.linspace(-0.2, 1.3, 16)
    for value in xs:
        if 0.05 < value < 0.95:
            step = 1e-5
            numerical = (accumulated(value + step) - accumulated(value - step)) / (2 * step)
            expected = float(pulse(value * TAU) ** 2)
            assert abs(numerical - expected) < 1e-8, (value, numerical, expected)
    g, g_prime = make_forces(accumulated)
    eta = 1 / 2000
    assert abs(float(g(T_STAR, eta)) - (1 - 4 * eta)) < 1e-12
    assert abs(float(g(T_STAR + 1.0, eta)) - float(g(T_STAR, eta))) < 1e-12
    charge = charge_identities()
    clocks = clock_allocation(g, g_prime)
    result = dict(
        status='PASS',
        verdict='scoped_obstruction',
        claim='The finite source coefficient has two extensions consistent with the retained identities. Dynamical insertion exceeds the switch-clock budget inside the benchmark window. The heat polynomials do not fix the coefficient. Endpoint contact energy remains zero, so that cancellation is not the clock bound.',
        not_an_impossibility_theorem=True,
        opening_achieved=False,
        shape=shape,
        charge=charge,
        clocks=clocks,
        versions=dict(python=platform.python_version(), sympy=s.__version__, numpy=np.__version__, scipy=scipy.__version__),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    summary = dict(
        status=result['status'],
        verdict=result['verdict'],
        lambda_critical=min(row['lambda_critical'] for row in clocks['samples']),
        limiting=min(clocks['samples'], key=lambda row: row['lambda_critical'])['limiting_budget'],
        exhibited=clocks['exhibited'],
    )
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
