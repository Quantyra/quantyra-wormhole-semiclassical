"""Admissible range of the dynamical source-clock coupling on the benchmark sample."""
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
ETA_MAX = 1 / 2000


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def decision_identities():
    """Opening clock varies y. A prescribed angle then contributes its charge."""
    y, A, rate, Nfield = s.symbols('y A rate N', positive=True)
    zeta, ell = s.symbols('zeta_s ell', positive=True)
    contact = A * y**2 * rate**2
    hamiltonian = s.simplify(y * s.diff(contact, y) - contact)
    charge = s.simplify(hamiltonian / y)
    zero(hamiltonian - contact)
    zero(charge - A * y * rate**2)
    lam = s.symbols('lambda', positive=True)
    zero(s.simplify(lam - A / (2 * Nfield)).subs(A, 2 * s.pi * zeta * ell**2).subs(lam, s.pi * zeta * ell**2 / Nfield))
    # The linear source carries eta**2. Smaller eta is not automatically safer once z changes,
    # so the checker compares two values rather than assuming the scaling.
    eta = s.symbols('eta', positive=True)
    zero(s.expand((4 * eta * rate) ** 2 - eta**2 * (4 * rate) ** 2))
    return dict(
        opening_rule='D',
        fixed_metric_statement='retained for prescribed theta and a metric that is not varied',
        varied_charge='A*y*(dTheta/dt)^2',
        lambda_definition='pi*zeta_s*ell**2/N',
        linear_source_scales_as='eta**2',
    )


def switch_grid():
    grid = np.linspace(0.0, alloc.T_STAR, 801)
    grid[-1] = alloc.T_STAR
    return grid


def evaluate(force, z0, v0, eta, lam, grid, base=None):
    if base is None:
        base = solve_ivp(
            alloc.clock_rhs(force, eta, 0.0), [0.0, alloc.T_END], [z0, v0],
            method='DOP853', rtol=1e-11, atol=1e-12, dense_output=True, t_eval=grid,
        )
        if not base.success:
            raise RuntimeError(base.message)
    finite = solve_ivp(
        alloc.clock_rhs(force, eta, lam), [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, dense_output=True, t_eval=grid,
    )
    if not finite.success:
        raise RuntimeError(finite.message)
    delta = finite.y - base.y
    index = -1
    gf = 1 - 4 * eta
    e_finite = float(finite.y[1, index] ** 2 + finite.y[0, index] ** 2 - gf * finite.y[0, index])
    e_base = float(base.y[1, index] ** 2 + base.y[0, index] ** 2 - gf * base.y[0, index])
    measured = dict(
        sup_abs_delta_z_on_switch=float(np.max(np.abs(delta[0]))),
        post_switch_amplitude=float(np.hypot(delta[0, index], delta[1, index])),
        abs_delta_e_f=abs(e_finite - e_base),
    )
    measured['admissible'] = bool(
        measured['sup_abs_delta_z_on_switch'] < SWITCH
        and measured['post_switch_amplitude'] < AFTER
        and measured['abs_delta_e_f'] < ENERGY
    )
    return measured, base


def linear_threshold(force, z0, v0, eta):
    grid = switch_grid()
    free = solve_ivp(
        alloc.clock_rhs(force, eta, 0.0), [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not free.success:
        raise RuntimeError(free.message)
    linear = alloc.integrate_linear(free.y[0], grid, eta, 'DOP853')
    index = -1
    switch = grid <= alloc.T_STAR + 1e-12
    sup_switch = float(np.max(np.abs(linear[0, switch])))
    amplitude = float(np.hypot(linear[0, index], linear[1, index]))
    gf = 1 - 4 * eta
    de = float(
        2 * free.y[1, index] * linear[1, index]
        + 2 * free.y[0, index] * linear[0, index]
        - gf * linear[0, index]
    )
    limits = [SWITCH / sup_switch, AFTER / amplitude, ENERGY / abs(de)]
    return float(min(limits)), sup_switch


def bracket(force, z0, v0, eta, sign, base, grid):
    """Largest |lambda| of this sign whose three clock budgets still hold."""
    lo = 0.0
    hi = 2000.0
    while True:
        measured, _ = evaluate(force, z0, v0, eta, sign * hi, grid, base)
        if not measured['admissible']:
            break
        hi *= 2.0
        if hi > 1.0e7:
            raise RuntimeError('budget remained admissible past 1e7')
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        measured, _ = evaluate(force, z0, v0, eta, sign * mid, grid, base)
        if measured['admissible']:
            lo = mid
        else:
            hi = mid
    edge, _ = evaluate(force, z0, v0, eta, sign * lo, grid, base)
    return lo, edge


def window(force):
    grid = switch_grid()
    corners = []
    for z0 in (0.25, 7 / 25):
        for v0 in (-ETA_MAX, 0.0, ETA_MAX):
            linear, response = linear_threshold(force, z0, v0, ETA_MAX)
            base_measured, base = evaluate(force, z0, v0, ETA_MAX, 0.0, grid)
            assert base_measured['admissible']
            positive, positive_edge = bracket(force, z0, v0, ETA_MAX, 1.0, base, grid)
            negative, negative_edge = bracket(force, z0, v0, ETA_MAX, -1.0, base, grid)
            half_linear, _ = linear_threshold(force, z0, v0, ETA_MAX / 2)
            corners.append(dict(
                eta=ETA_MAX, z0=z0, v0=v0,
                linear_lambda_threshold=linear,
                linear_switch_response=response,
                nonlinear_positive_limit=positive,
                nonlinear_negative_limit=-negative,
                positive_edge=positive_edge,
                negative_edge=negative_edge,
                half_eta_linear_lambda_threshold=half_linear,
            ))
    positive_limit = min(row['nonlinear_positive_limit'] for row in corners)
    negative_limit = max(row['nonlinear_negative_limit'] for row in corners)
    worst = min(corners, key=lambda row: row['nonlinear_positive_limit'])
    # The half-eta linear threshold must be larger at every sampled corner.
    for row in corners:
        assert row['half_eta_linear_lambda_threshold'] > row['linear_lambda_threshold']
    # Confirm the worst positive edge with an independent integrator.
    other = solve_ivp(
        alloc.clock_rhs(force, ETA_MAX, worst['nonlinear_positive_limit']),
        [0.0, alloc.T_END], [worst['z0'], worst['v0']],
        method='RK45', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not other.success:
        raise RuntimeError(other.message)
    dop = solve_ivp(
        alloc.clock_rhs(force, ETA_MAX, worst['nonlinear_positive_limit']),
        [0.0, alloc.T_END], [worst['z0'], worst['v0']],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    gap = float(np.max(np.abs(other.y - dop.y)))
    assert gap < 1e-8, gap
    # The bisection assumes the switch deviation grows with |lambda|.
    # Check that growth on the worst corner before calling the interval admissible.
    _, worst_base = evaluate(force, worst['z0'], worst['v0'], ETA_MAX, 0.0, grid)
    previous = 0.0
    for sign, limit in ((1.0, worst['nonlinear_positive_limit']), (-1.0, abs(worst['nonlinear_negative_limit']))):
        for fraction in (0.25, 0.5, 0.75, 1.0):
            probed, _ = evaluate(force, worst['z0'], worst['v0'], ETA_MAX, sign * fraction * limit, grid, worst_base)
            assert probed['admissible'], (sign, fraction, probed)
            assert probed['sup_abs_delta_z_on_switch'] >= previous - 1e-8
            previous = probed['sup_abs_delta_z_on_switch']
        previous = 0.0
    plain = []
    for row in corners:
        plain.append(dict(
            eta=row['eta'], z0=row['z0'], v0=row['v0'],
            linear_lambda_threshold=row['linear_lambda_threshold'],
            nonlinear_positive_limit=row['nonlinear_positive_limit'],
            nonlinear_negative_limit=row['nonlinear_negative_limit'],
            half_eta_linear_lambda_threshold=row['half_eta_linear_lambda_threshold'],
            positive_edge_switch=row['positive_edge']['sup_abs_delta_z_on_switch'],
            positive_edge_after=row['positive_edge']['post_switch_amplitude'],
            positive_edge_energy=row['positive_edge']['abs_delta_e_f'],
        ))
    return dict(
        sample='eta=1/2000 and the six initial corners of the retained benchmark',
        lambda_interval=[negative_limit, positive_limit],
        zeta_s_interval='lambda*N/(pi*ell**2)',
        no_retained_input_selects_a_point=True,
        area_budget='not evaluated; Phi still depends on the unallocated y and E',
        worst_positive_corner=dict(z0=worst['z0'], v0=worst['v0'], lambda_limit=worst['nonlinear_positive_limit']),
        independent_integrator_gap=gap,
        corners=plain,
    )


def main():
    warnings.filterwarnings('error', category=scipy.integrate.IntegrationWarning)
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    primitive = alloc.shape_identities().pop('primitive_polynomial')
    force, _ = alloc.make_forces(alloc.antiderivative_factory(primitive))
    result = dict(
        status='PASS',
        verdict='admissible_range_on_benchmark_sample',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        decision=decision_identities(),
        window=window(force),
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
        rule=result['decision']['opening_rule'],
        lambda_interval=result['window']['lambda_interval'],
        worst=result['window']['worst_positive_corner'],
        s074_gate_closed=False,
    )
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
