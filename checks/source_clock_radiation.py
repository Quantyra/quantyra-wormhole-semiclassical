"""Identify Lhat on the varied clock and evaluate the radiation law on a realized angle."""
import argparse
import hashlib
import json
import platform
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import sympy as s
from scipy.integrate import solve_ivp
import scipy

import source_clock_allocation as alloc


SWITCH = alloc.SWITCH_BUDGET
AFTER = alloc.AFTER_BUDGET
ENERGY = alloc.ENERGY_BUDGET
ETA = 1 / 2000
POSITIVE_ANCHOR = 4875.4563331604
NEGATIVE_ANCHOR = -4757.61079788208
FRACTIONS = (0.25, 0.5, 0.75, 1.0)
OVERSHOOT = 1.1
OVERSHOOT_KAPPA = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
WORST_Z0 = 7 / 25
WORST_V0 = ETA
LAMBDAS = (
    0.0,
    0.5 * POSITIVE_ANCHOR,
    POSITIVE_ANCHOR,
    0.5 * NEGATIVE_ANCHOR,
    NEGATIVE_ANCHOR,
)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def geometry_identities():
    """Product-throat Lhat weight in the cutoff coordinate."""
    re, eps, sch, rate = s.symbols('r_e epsilon Sch rate', positive=True)
    pi = s.pi
    curvature = (1 / re) * (1 - eps**2 * sch)
    gamma_tt = -(re**2) / eps**2
    l_tt = curvature * gamma_tt
    trace = s.simplify(l_tt / gamma_tt)
    zero(trace - curvature)
    lhat_tt = s.simplify(l_tt - trace * gamma_tt / 3)
    zero(lhat_tt - s.Rational(2, 3) * curvature * gamma_tt)
    # Angular extrinsic curvature is zero, so the surface is not umbilic.
    lhat_sphere = s.simplify(0 - trace * re**2 / 3)
    assert s.simplify(lhat_sphere) != 0
    raised = s.simplify(lhat_tt / gamma_tt**2)
    measure = 4 * pi * re**3 / eps
    density = s.simplify(measure * raised * rate**2)
    expected = -s.Rational(8, 3) * pi * eps * (1 - eps**2 * sch) * rate**2
    zero(density - expected)
    series = density.series(eps, 0, 3).removeO()
    zero(series + s.Rational(8, 3) * pi * eps * rate**2)
    zero(density.subs(eps, 0))
    return dict(
        integrated_density='-(8*pi/3)*epsilon*(1-epsilon**2*Sch)*rate**2',
        conformal_boundary_weight=0,
        umbilic=False,
        schwarzian_order='epsilon**3',
        identified_input='O(epsilon) product-throat weight; absent from the JT clock',
    )


def delay_identity():
    theta, Q, eta = s.symbols('theta Q eta', real=True)
    energy = (-theta + Q) / (4 * s.pi)
    drive = -energy / eta
    zero(drive - (theta - Q) / (4 * s.pi * eta))
    # Returned-switch sum. Shift by one period and drop the image that has
    # not started; the difference is the current junction rate.
    rho0, rho1, rho2 = s.symbols('rho0 rho1 rho2', real=True)
    current = rho0 + rho1 + rho2
    previous = rho1 + rho2
    zero(current - previous - rho0)
    return dict(
        drive='(theta-Q)/(4*pi*eta)',
        delay='D(t)=rho(t)+D(t-pi)',
        prescribed_rate='4*eta*(S(t/tau)-S((t-pi)/tau))',
    )


def _delayed_D(t, history, eta, kappa):
    if t <= 1e-12:
        return 0.0
    for t0, t1, sol in history:
        if t0 <= t <= t1 + 1e-9:
            z = float(sol(min(t, t1))[0])
            r0 = 4.0 * eta * float(alloc.angle_rate_factor(t))
            return r0 + kappa * z + _delayed_D(t - np.pi, history, eta, kappa)
    raise RuntimeError('delay time outside recorded history: %s' % t)


def march(eta, lam, kappa, z0, v0, method):
    """Method of steps for the realized returned-switch drive."""
    history = []
    state = np.array([z0, v0, 0.0, 0.0], dtype=float)
    t0 = 0.0
    while t0 < alloc.T_END - 1e-12:
        t1 = min(t0 + np.pi, alloc.T_END)

        def rhs(t, y, history=history, eta=eta, lam=lam, kappa=kappa):
            z, zt, theta, charge = y
            r0 = 4.0 * eta * float(alloc.angle_rate_factor(t))
            rho = r0 + kappa * z
            previous = _delayed_D(t - np.pi, history, eta, kappa)
            profile = rho + previous
            drive = (theta - charge) / (4.0 * np.pi * eta)
            accel = (
                -(1.0 + lam * r0 * r0) * z
                - 4.0 * lam * kappa * r0 * z * z
                - 3.0 * lam * kappa * kappa * z * z * z
                + 0.5 * drive
            )
            return [zt, accel, rho, profile * profile - previous * previous]

        solution = solve_ivp(
            rhs, [t0, t1], state, method=method, rtol=1e-11, atol=1e-12, dense_output=True,
        )
        if not solution.success or not np.all(np.isfinite(solution.y)):
            return None
        history.append((t0, t1, solution.sol))
        state = solution.y[:, -1]
        t0 = t1
    return history


def sample_state(history, grid):
    out = np.empty((4, len(grid)), dtype=float)
    for i, time in enumerate(grid):
        found = False
        for t0, t1, sol in history:
            if t0 <= time <= t1 + 1e-9:
                out[:, i] = sol(min(time, t1))
                found = True
                break
        if not found:
            raise RuntimeError('sample time outside history')
    return out


def grids():
    switch = np.linspace(0.0, alloc.T_STAR, 801)
    switch[-1] = alloc.T_STAR
    post = np.linspace(alloc.T_STAR, alloc.T_END, 1601)
    post[0] = alloc.T_STAR
    return np.unique(np.concatenate([switch, post]))


def drive_regression(force):
    history = march(ETA, 0.0, 0.0, WORST_Z0, WORST_V0, 'DOP853')
    if history is None:
        raise RuntimeError('prescribed march failed')
    times = np.linspace(0.0, alloc.T_END, 401)
    times[-1] = alloc.T_END
    state = sample_state(history, times)
    errors = []
    for i, time in enumerate(times):
        drive = (state[2, i] - state[3, i]) / (4.0 * np.pi * ETA)
        target = float(np.asarray(force(time, ETA), dtype=float))
        errors.append(abs(drive - target))
    gap = float(max(errors))
    assert gap < 1e-8, gap
    star = sample_state(history, np.array([alloc.T_STAR]))
    final = (star[2, 0] - star[3, 0]) / (4.0 * np.pi * ETA)
    assert abs(final - (1.0 - 4.0 * ETA)) < 1e-8, final
    return dict(max_abs_drive_error=gap, final_drive=float(final), final_target=1.0 - 4.0 * ETA)


def anchor_regression(force):
    grid = np.linspace(0.0, alloc.T_END, 401)
    grid[-1] = alloc.T_END
    history = march(ETA, POSITIVE_ANCHOR, 0.0, WORST_Z0, WORST_V0, 'DOP853')
    if history is None:
        raise RuntimeError('anchor march failed')
    mine = sample_state(history, grid)[:2]
    retained = solve_ivp(
        alloc.clock_rhs(force, ETA, POSITIVE_ANCHOR),
        [0.0, alloc.T_END], [WORST_Z0, WORST_V0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not retained.success:
        raise RuntimeError(retained.message)
    gap = float(np.max(np.abs(mine - retained.y)))
    assert gap < 1e-6, gap
    return dict(anchor_trajectory_gap=gap, lambda_value=POSITIVE_ANCHOR)


def corners():
    return [(z0, v0) for z0 in (0.25, WORST_Z0) for v0 in (-ETA, 0.0, ETA)]


def measure(history_or_none, base, grid, eta):
    if history_or_none is None:
        return dict(admissible=False, left_domain=True, integrator_failed=True)
    finite = sample_state(history_or_none, grid)
    if not np.all(np.isfinite(finite[:2])):
        return dict(admissible=False, left_domain=True, integrator_failed=True)
    star = int(np.searchsorted(grid, alloc.T_STAR))
    delta = finite[:2] - base
    gf = 1.0 - 4.0 * eta
    e_finite = float(finite[1, star] ** 2 + finite[0, star] ** 2 - gf * finite[0, star])
    e_base = float(base[1, star] ** 2 + base[0, star] ** 2 - gf * base[0, star])
    post = grid > alloc.T_STAR
    measured = dict(
        sup_abs_delta_z_on_switch=float(np.max(np.abs(delta[0, :star + 1]))),
        post_switch_sup_abs_delta_z=float(np.max(np.abs(delta[0, post]))),
        abs_delta_e_f=abs(e_finite - e_base),
        min_z=float(np.min(finite[0])),
        left_domain=bool(np.min(finite[0]) <= 0.0),
        integrator_failed=False,
    )
    measured['admissible'] = bool(
        measured['sup_abs_delta_z_on_switch'] < SWITCH
        and measured['post_switch_sup_abs_delta_z'] < AFTER
        and measured['abs_delta_e_f'] < ENERGY
        and not measured['left_domain']
    )
    ratios = {
        'switch': measured['sup_abs_delta_z_on_switch'] / SWITCH,
        'post_switch': measured['post_switch_sup_abs_delta_z'] / AFTER,
        'energy': measured['abs_delta_e_f'] / ENERGY,
    }
    measured['binding_budget'] = max(ratios, key=ratios.get)
    return measured


def bracket_kappa(z0, v0, eta, lam, sign, grid, base):
    lo = 0.0
    hi = 1e-4
    while True:
        measured = measure(march(eta, lam, sign * hi, z0, v0, 'DOP853'), base, grid, eta)
        if not measured['admissible']:
            break
        hi *= 2.0
        if hi > 1.0e3:
            raise RuntimeError('radiation budget remained admissible past 1e3')
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        measured = measure(march(eta, lam, sign * mid, z0, v0, 'DOP853'), base, grid, eta)
        if measured['admissible']:
            lo = mid
        else:
            hi = mid
    edge = measure(march(eta, lam, sign * lo, z0, v0, 'DOP853'), base, grid, eta)
    return lo, edge


def window():
    grid = grids()
    bases = {}
    for z0, v0 in corners():
        history = march(ETA, 0.0, 0.0, z0, v0, 'DOP853')
        if history is None:
            raise RuntimeError('base march failed')
        bases[(z0, v0)] = sample_state(history, grid)[:2]
    samples = []
    for lam in LAMBDAS:
        positive = []
        negative = []
        for z0, v0 in corners():
            base = bases[(z0, v0)]
            pos, pos_edge = bracket_kappa(z0, v0, ETA, lam, 1.0, grid, base)
            neg, neg_edge = bracket_kappa(z0, v0, ETA, lam, -1.0, grid, base)
            positive.append(dict(z0=z0, v0=v0, kappa=pos, edge=pos_edge))
            negative.append(dict(z0=z0, v0=v0, kappa=-neg, edge=neg_edge))
            print('corner', lam, z0, v0, -neg, pos, flush=True)
        shared_positive = min(row['kappa'] for row in positive)
        shared_negative = max(row['kappa'] for row in negative)
        pos_worst = min(positive, key=lambda row: row['kappa'])
        neg_worst = max(negative, key=lambda row: row['kappa'])
        print('sample', lam, shared_negative, shared_positive, pos_worst['edge'].get('binding_budget'), flush=True)
        samples.append(dict(
            lambda_value=lam,
            kappa_positive=shared_positive,
            kappa_negative=shared_negative,
            positive_binding_budget=pos_worst['edge'].get('binding_budget'),
            negative_binding_budget=neg_worst['edge'].get('binding_budget'),
            positive_edge=dict(
                z0=pos_worst['z0'], v0=pos_worst['v0'],
                switch=pos_worst['edge'].get('sup_abs_delta_z_on_switch'),
                post_switch=pos_worst['edge'].get('post_switch_sup_abs_delta_z'),
                energy=pos_worst['edge'].get('abs_delta_e_f'),
                admissible=pos_worst['edge'].get('admissible'),
            ),
            negative_edge=dict(
                z0=neg_worst['z0'], v0=neg_worst['v0'],
                switch=neg_worst['edge'].get('sup_abs_delta_z_on_switch'),
                post_switch=neg_worst['edge'].get('post_switch_sup_abs_delta_z'),
                energy=neg_worst['edge'].get('abs_delta_e_f'),
                admissible=neg_worst['edge'].get('admissible'),
            ),
        ))
    return samples


def monotonicity(samples):
    sample = next(row for row in samples if row['lambda_value'] == 0.0)
    grid = grids()
    probed = []
    endpoints = (
        (1.0, sample['kappa_positive'], sample['positive_edge']),
        (-1.0, abs(sample['kappa_negative']), sample['negative_edge']),
    )
    for sign, limit, edge in endpoints:
        z0 = edge['z0']
        v0 = edge['v0']
        base = sample_state(march(ETA, 0.0, 0.0, z0, v0, 'DOP853'), grid)[:2]
        for fraction in FRACTIONS:
            measured = measure(
                march(ETA, 0.0, sign * fraction * limit, z0, v0, 'DOP853'),
                base, grid, ETA,
            )
            assert measured['admissible'], (sign, fraction, z0, v0, measured)
            probed.append(dict(
                sign=sign, fraction=fraction, z0=z0, v0=v0, admissible=True,
                switch=measured['sup_abs_delta_z_on_switch'],
                post_switch=measured['post_switch_sup_abs_delta_z'],
                energy=measured['abs_delta_e_f'],
            ))
    return dict(lambda_value=0.0, probe=probed)


def independent_edge(samples):
    sample = next(row for row in samples if row['lambda_value'] == 0.0)
    grid = grids()
    z0 = sample['positive_edge']['z0']
    v0 = sample['positive_edge']['v0']
    kappa = sample['kappa_positive']
    dop = sample_state(march(ETA, 0.0, kappa, z0, v0, 'DOP853'), grid)[:2]
    other = sample_state(march(ETA, 0.0, kappa, z0, v0, 'RK45'), grid)[:2]
    gap = float(np.max(np.abs(dop - other)))
    assert gap < 1e-7, gap
    return dict(lambda_value=0.0, kappa=kappa, corner=dict(z0=z0, v0=v0), dop853_rk45_gap=gap)


def overshoot_probe():
    grid = grids()
    base = sample_state(march(ETA, 0.0, 0.0, WORST_Z0, WORST_V0, 'DOP853'), grid)[:2]
    lam = OVERSHOOT * POSITIVE_ANCHOR
    rows = []
    zero_kappa = measure(march(ETA, lam, 0.0, WORST_Z0, WORST_V0, 'DOP853'), base, grid, ETA)
    for magnitude in OVERSHOOT_KAPPA:
        for sign in (1.0, -1.0):
            measured = measure(
                march(ETA, lam, sign * magnitude, WORST_Z0, WORST_V0, 'DOP853'),
                base, grid, ETA,
            )
            rows.append(dict(
                kappa=sign * magnitude,
                admissible=measured['admissible'],
                binding_budget=measured.get('binding_budget'),
                sup_abs_delta_z_on_switch=measured.get('sup_abs_delta_z_on_switch'),
                post_switch_sup_abs_delta_z=measured.get('post_switch_sup_abs_delta_z'),
                abs_delta_e_f=measured.get('abs_delta_e_f'),
            ))
    return dict(
        lambda_overshoot=lam,
        corner=dict(z0=WORST_Z0, v0=WORST_V0),
        kappa_zero_admissible=zero_kappa['admissible'],
        samples=rows,
        any_sample_admissible=any(row['admissible'] for row in rows),
    )


def main():
    warnings.filterwarnings('error', category=scipy.integrate.IntegrationWarning)
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    primitive = alloc.shape_identities().pop('primitive_polynomial')
    force, _ = alloc.make_forces(alloc.antiderivative_factory(primitive))
    geometry = geometry_identities()
    delay = delay_identity()
    drive = drive_regression(force)
    anchor = anchor_regression(force)
    samples = window()
    mono = monotonicity(samples)
    edge = independent_edge(samples)
    overshoot = overshoot_probe()
    result = dict(
        status='PASS',
        verdict='lhat_weight_vanishes_and_radiation_kappa_range',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        area_budget='not evaluated; no replacement target was written before the test',
        geometry=geometry,
        delay=delay,
        drive_regression=drive,
        anchor_regression=anchor,
        c_theta_separate_range=False,
        radiation_ranges=samples,
        monotonicity=mono,
        independent_integrator=edge,
        overshoot_probe=overshoot,
        versions=dict(
            python=platform.python_version(), sympy=s.__version__,
            numpy=np.__version__, scipy=scipy.__version__,
        ),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    summary = dict(
        status=result['status'],
        verdict=result['verdict'],
        drive_error=drive['max_abs_drive_error'],
        anchor_gap=anchor['anchor_trajectory_gap'],
        kappa_at_zero=[samples[0]['kappa_negative'], samples[0]['kappa_positive']],
        s074_gate_closed=False,
        opening_achieved=False,
    )
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
