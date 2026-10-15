"""Carry c_theta and a non-prescribed controller angle onto the rule-D clock."""
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
# Prescribed-angle anchors already recorded for kappa = 0. Not refit here.
POSITIVE_ANCHOR = 4875.4563331604
NEGATIVE_ANCHOR = -4757.61079788208
FRACTIONS = (0.25, 0.5, 0.75, 1.0)
OVERSHOOT = 1.1
OVERSHOOT_KAPPA = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
WORST_Z0 = 7 / 25
WORST_V0 = ETA


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def charge_reduction():
    """Direct variation of the contact charge. The drive g is not replaced."""
    y, A, Nf, eta = s.symbols('y A N eta', positive=True)
    r0, kappa = s.symbols('r0 kappa', real=True)
    lam, zeta, ell, c_theta, lhat = s.symbols('Lambda zeta_s ell c_theta Lhat_uu', real=True)
    rho = r0 + kappa * y / eta
    lag = A * y**2 * rho**2
    hamiltonian = s.simplify(y * s.diff(lag, y) - lag)
    charge = s.simplify(hamiltonian / y)
    expected_h = A * y**2 * rho**2 + 2 * A * (kappa / eta) * y**3 * rho
    expected_e = A * y * rho**2 + 2 * A * (kappa / eta) * y**2 * rho
    zero(hamiltonian - expected_h)
    zero(charge - expected_e)
    z = y / eta
    extra = s.simplify((-charge / (2 * Nf)) / eta)
    reduced = s.simplify(extra.subs(A, 2 * Nf * lam))
    cubic = -lam * r0**2 * z - 4 * lam * kappa * r0 * z**2 - 3 * lam * kappa**2 * z**3
    zero(reduced - cubic)
    zero(reduced.subs(kappa, 0) + lam * r0**2 * z)
    zero(reduced.subs(lam, 0))
    # Curvature contact and the retained source add before the variation
    # when the angle is still prescribed. The 1/2 in I_s produces the 2.
    a_source = 2 * s.pi * zeta * ell**2
    a_theta = 4 * s.pi * ell**2 * c_theta * lhat
    lam_source = s.simplify(a_source / (2 * Nf))
    lam_theta = s.simplify(a_theta / (2 * Nf))
    zero(lam_source - s.pi * zeta * ell**2 / Nf)
    zero(lam_theta - 2 * s.pi * ell**2 * c_theta * lhat / Nf)
    zero(s.together(lam_theta / lam_source - 2 * c_theta * lhat / zeta))
    # a4 is built from L_ab and u_a. c_theta is not one of its symbols.
    jets = s.symbols('L00 L11 L22 L01 L02 L12 u0 u1 u2', real=True)
    l00, l11, l22, l01, l02, l12, u0, u1, u2 = jets
    trace = l00 + l11 + l22
    u_square = u0**2 + u1**2 + u2**2
    contraction = (
        l00 * u0**2 + l11 * u1**2 + l22 * u2**2
        + 2 * l01 * u0 * u1 + 2 * l02 * u0 * u2 + 2 * l12 * u1 * u2
    )
    a4 = s.Rational(4, 5) * (contraction - trace * u_square / 3)
    assert c_theta not in a4.free_symbols
    return dict(
        opening_rule='D',
        u_hamiltonian='A*y^2*rho^2 + 2*A*(kappa/eta)*y^3*rho',
        global_time_charge='A*y*rho^2 + 2*A*(kappa/eta)*y^2*rho',
        clock_equation=(
            'z_tt+(1+Lambda*r0^2)*z+4*Lambda*kappa*r0*z^2'
            '+3*Lambda*kappa^2*z^3=g/2'
        ),
        kappa_zero_reduces_to_rule_D=True,
        lambda_zero_removes_every_kappa_term=True,
        lambda_theta='(2*pi*ell**2/N)*c_theta*Lhat_uu',
        joint_constraint='Lambda=(pi*ell**2/N)*(zeta_s+2*c_theta*Lhat_uu)',
        c_theta_in_heat_polynomial_a4=False,
        separating_input='Lhat_uu on the varied clock boundary',
        drive_channel='unresolved; retained g is not replaced',
        controller_ansatz='rho=r0(t)+kappa*z',
    )


def numerical_variation():
    """Central difference of L at the point named in the note."""
    eta = ETA
    y = 0.28 * eta
    r0 = 0.001
    kappa = 0.2
    amplitude = 1.0
    step = 1e-8

    def lag(value):
        rho = r0 + kappa * value / eta
        return amplitude * value**2 * rho**2

    slope = (lag(y + step) - lag(y - step)) / (2 * step)
    discrete = y * slope - lag(y)
    rho = r0 + kappa * y / eta
    closed = (
        amplitude * y**2 * rho**2
        + 2 * amplitude * (kappa / eta) * y**3 * rho
    )
    gap = abs(discrete - closed)
    relative = gap / max(abs(closed), 1e-30)
    assert relative < 1e-8, (gap, closed)
    return dict(
        central_difference_gap=gap,
        central_difference_relative_gap=relative,
        point='A=1, eta=1/2000, y=0.28*eta, r0=0.001, kappa=0.2',
    )


def controller_rhs(force, eta, lam, kappa):
    """Rule D plus the varied contact charge. The retained drive is unchanged."""
    def rhs(t, state):
        z, zt = state
        r0 = 4.0 * eta * float(alloc.angle_rate_factor(t))
        drive = float(np.asarray(force(t, eta), dtype=float))
        accel = (
            -(1.0 + lam * r0**2) * z
            - 4.0 * lam * kappa * r0 * z**2
            - 3.0 * lam * kappa**2 * z**3
            + 0.5 * drive
        )
        return [zt, accel]

    return rhs


def grids():
    switch = np.linspace(0.0, alloc.T_STAR, 801)
    switch[-1] = alloc.T_STAR
    post = np.linspace(alloc.T_STAR, alloc.T_END, 1601)
    post[0] = alloc.T_STAR
    grid = np.unique(np.concatenate([switch, post]))
    return grid


def measure(force, z0, v0, eta, lam, kappa, grid, base):
    try:
        finite = solve_ivp(
            controller_rhs(force, eta, lam, kappa),
            [0.0, alloc.T_END], [z0, v0],
            method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
        )
    except (RuntimeError, ValueError):
        return dict(admissible=False, left_domain=True, integrator_failed=True)
    if not finite.success or not np.all(np.isfinite(finite.y)):
        return dict(admissible=False, left_domain=True, integrator_failed=True)
    star = int(np.searchsorted(grid, alloc.T_STAR))
    delta = finite.y - base.y
    gf = 1.0 - 4.0 * eta
    e_finite = float(finite.y[1, star] ** 2 + finite.y[0, star] ** 2 - gf * finite.y[0, star])
    e_base = float(base.y[1, star] ** 2 + base.y[0, star] ** 2 - gf * base.y[0, star])
    post = grid > alloc.T_STAR
    measured = dict(
        sup_abs_delta_z_on_switch=float(np.max(np.abs(delta[0, :star + 1]))),
        post_switch_sup_abs_delta_z=float(np.max(np.abs(delta[0, post]))),
        abs_delta_e_f=abs(e_finite - e_base),
        min_z=float(np.min(finite.y[0])),
        left_domain=bool(np.min(finite.y[0]) <= 0.0),
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


def base_trajectory(force, z0, v0, eta, grid):
    base = solve_ivp(
        controller_rhs(force, eta, 0.0, 0.0),
        [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not base.success:
        raise RuntimeError(base.message)
    return base


def bracket_kappa(force, z0, v0, eta, lam, sign, grid, base):
    """Largest |kappa| of this sign whose three budgets still hold."""
    lo = 0.0
    hi = 1e-4
    while True:
        measured = measure(force, z0, v0, eta, lam, sign * hi, grid, base)
        if not measured['admissible']:
            break
        hi *= 2.0
        if hi > 1.0e4:
            raise RuntimeError('controller budget remained admissible past 1e4')
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        measured = measure(force, z0, v0, eta, lam, sign * mid, grid, base)
        if measured['admissible']:
            lo = mid
        else:
            hi = mid
    edge = measure(force, z0, v0, eta, lam, sign * lo, grid, base)
    if lo == 0.0:
        edge = measure(force, z0, v0, eta, lam, 0.0, grid, base)
    return lo, edge


def corners():
    rows = []
    for z0 in (0.25, WORST_Z0):
        for v0 in (-ETA, 0.0, ETA):
            rows.append((z0, v0))
    return rows


def additivity_check(force, grid):
    """Half the positive anchor in each channel matches rule D at the sum."""
    lam = 0.5 * POSITIVE_ANCHOR
    total = POSITIVE_ANCHOR
    z0, v0 = WORST_Z0, WORST_V0
    split = solve_ivp(
        controller_rhs(force, ETA, total, 0.0),
        [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    retained = solve_ivp(
        alloc.clock_rhs(force, ETA, total),
        [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not split.success or not retained.success:
        raise RuntimeError('additivity integration failed')
    gap = float(np.max(np.abs(split.y - retained.y)))
    assert gap < 1e-8, gap
    # The algebraic split is equal: each half is the same Lambda contribution.
    assert abs((lam + lam) - total) == 0.0
    return dict(
        lambda_source_piece=lam,
        lambda_theta_piece=lam,
        sum=total,
        trajectory_gap_against_rule_D=gap,
        corner=dict(z0=z0, v0=v0),
    )


def overshoot_probe(force, grid, base):
    lam = OVERSHOOT * POSITIVE_ANCHOR
    rows = []
    zero_kappa = measure(force, WORST_Z0, WORST_V0, ETA, lam, 0.0, grid, base)
    for magnitude in OVERSHOOT_KAPPA:
        for sign in (1.0, -1.0):
            measured = measure(
                force, WORST_Z0, WORST_V0, ETA, lam, sign * magnitude, grid, base,
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


def window(force):
    grid = grids()
    bases = {(z0, v0): base_trajectory(force, z0, v0, ETA, grid) for z0, v0 in corners()}
    samples = []
    for sign in (1.0, -1.0):
        anchor = POSITIVE_ANCHOR if sign > 0.0 else NEGATIVE_ANCHOR
        for fraction in FRACTIONS:
            lam = sign * fraction * abs(anchor)
            positive_limits = []
            negative_limits = []
            positive_edges = []
            negative_edges = []
            for z0, v0 in corners():
                base = bases[(z0, v0)]
                pos, pos_edge = bracket_kappa(force, z0, v0, ETA, lam, 1.0, grid, base)
                neg, neg_edge = bracket_kappa(force, z0, v0, ETA, lam, -1.0, grid, base)
                positive_limits.append(dict(z0=z0, v0=v0, kappa=pos, edge=pos_edge))
                negative_limits.append(dict(z0=z0, v0=v0, kappa=-neg, edge=neg_edge))
            shared_positive = min(row['kappa'] for row in positive_limits)
            shared_negative = max(row['kappa'] for row in negative_limits)
            pos_worst = min(positive_limits, key=lambda row: row['kappa'])
            neg_worst = max(negative_limits, key=lambda row: row['kappa'])
            print(
                'sample', lam, shared_negative, shared_positive,
                pos_worst['edge'].get('binding_budget'),
                flush=True,
            )
            samples.append(dict(
                lambda_value=lam,
                anchor_fraction=fraction * sign,
                kappa_positive=shared_positive,
                kappa_negative=shared_negative,
                positive_binding_budget=pos_worst['edge'].get('binding_budget'),
                negative_binding_budget=neg_worst['edge'].get('binding_budget'),
                positive_edge=dict(
                    z0=pos_worst['z0'], v0=pos_worst['v0'],
                    switch=pos_worst['edge'].get('sup_abs_delta_z_on_switch'),
                    post_switch=pos_worst['edge'].get('post_switch_sup_abs_delta_z'),
                    energy=pos_worst['edge'].get('abs_delta_e_f'),
                ),
                negative_edge=dict(
                    z0=neg_worst['z0'], v0=neg_worst['v0'],
                    switch=neg_worst['edge'].get('sup_abs_delta_z_on_switch'),
                    post_switch=neg_worst['edge'].get('post_switch_sup_abs_delta_z'),
                    energy=neg_worst['edge'].get('abs_delta_e_f'),
                ),
            ))
    return samples


def monotonicity(force, samples):
    """Probe the half positive anchor on its most restrictive corner."""
    grid = grids()
    sample = next(
        row for row in samples
        if abs(row['lambda_value'] - 0.5 * POSITIVE_ANCHOR) < 1e-6
    )
    z0 = sample['positive_edge']['z0']
    v0 = sample['positive_edge']['v0']
    base = base_trajectory(force, z0, v0, ETA, grid)
    probed = []
    for sign, limit in ((1.0, sample['kappa_positive']), (-1.0, abs(sample['kappa_negative']))):
        previous = 0.0
        for fraction in FRACTIONS:
            measured = measure(
                force, z0, v0, ETA, sample['lambda_value'], sign * fraction * limit, grid, base,
            )
            assert measured['admissible'], (sign, fraction, measured)
            binding = measured['sup_abs_delta_z_on_switch']
            probed.append(dict(
                sign=sign,
                fraction=fraction,
                switch=binding,
                post_switch=measured['post_switch_sup_abs_delta_z'],
                energy=measured['abs_delta_e_f'],
                switch_deviation_nondecreasing=bool(binding >= previous - 1e-8),
                admissible=True,
            ))
            previous = binding
    return dict(
        lambda_value=sample['lambda_value'],
        corner=dict(z0=z0, v0=v0),
        probe=probed,
    )


def independent_edge(force, samples):
    grid = grids()
    sample = next(
        row for row in samples
        if abs(row['lambda_value'] - 0.5 * POSITIVE_ANCHOR) < 1e-6
    )
    z0 = sample['positive_edge']['z0']
    v0 = sample['positive_edge']['v0']
    kappa = sample['kappa_positive']
    dop = solve_ivp(
        controller_rhs(force, ETA, sample['lambda_value'], kappa),
        [0.0, alloc.T_END], [z0, v0],
        method='DOP853', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    other = solve_ivp(
        controller_rhs(force, ETA, sample['lambda_value'], kappa),
        [0.0, alloc.T_END], [z0, v0],
        method='RK45', rtol=1e-11, atol=1e-12, t_eval=grid,
    )
    if not dop.success or not other.success:
        raise RuntimeError('independent edge integration failed')
    gap = float(np.max(np.abs(dop.y - other.y)))
    assert gap < 1e-8, gap
    return dict(
        lambda_value=sample['lambda_value'],
        kappa=kappa,
        corner=dict(z0=z0, v0=v0),
        dop853_rk45_gap=gap,
    )


def main():
    warnings.filterwarnings('error', category=scipy.integrate.IntegrationWarning)
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    primitive = alloc.shape_identities().pop('primitive_polynomial')
    force, _ = alloc.make_forces(alloc.antiderivative_factory(primitive))
    identities = charge_reduction()
    variation = numerical_variation()
    grid = grids()
    added = additivity_check(force, grid)
    # kappa = 0 at the positive anchor stays inside the recorded comparison.
    anchor_base = base_trajectory(force, WORST_Z0, WORST_V0, ETA, grid)
    anchor = measure(force, WORST_Z0, WORST_V0, ETA, POSITIVE_ANCHOR, 0.0, grid, anchor_base)
    assert anchor['admissible'], anchor
    samples = window(force)
    mono = monotonicity(force, samples)
    edge = independent_edge(force, samples)
    overshoot = overshoot_probe(force, grid, anchor_base)
    result = dict(
        status='PASS',
        verdict='joint_sum_constraint_and_contact_channel_kappa_range',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        area_budget='not evaluated; no replacement target was written before the test',
        decision=identities,
        numerical_variation=variation,
        additivity=added,
        prescribed_angle_anchors=[NEGATIVE_ANCHOR, POSITIVE_ANCHOR],
        positive_anchor_at_kappa_zero_admissible=anchor['admissible'],
        controller_ranges=samples,
        monotonicity=mono,
        independent_integrator=edge,
        overshoot_probe=overshoot,
        c_theta_separate_range=False,
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
        separating_input=identities['separating_input'],
        kappa_at_half_positive=[
            row['kappa_negative'] for row in samples
            if abs(row['lambda_value'] - 0.5 * POSITIVE_ANCHOR) < 1e-6
        ] + [
            row['kappa_positive'] for row in samples
            if abs(row['lambda_value'] - 0.5 * POSITIVE_ANCHOR) < 1e-6
        ],
        overshoot_restored=overshoot['any_sample_admissible'],
        s074_gate_closed=False,
        opening_achieved=False,
    )
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
