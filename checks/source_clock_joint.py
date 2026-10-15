"""Range the finite Weyl product on the realized returned-switch angle."""
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
import scipy

import source_clock_allocation as alloc
import source_clock_radiation as rad


RECORDED_KAPPA = {
    0.0: (-5.16183614730835e-05, 5.072014927864075e-05),
    0.5 * rad.POSITIVE_ANCHOR: (-4.602195620536803e-05, 5.619207024574278e-05),
    0.5 * rad.NEGATIVE_ANCHOR: (-5.705349445343019e-05, 4.523511528968811e-05),
}
PROBE_LAMBDAS = (
    0.5 * rad.POSITIVE_ANCHOR,
    rad.POSITIVE_ANCHOR,
    0.5 * rad.NEGATIVE_ANCHOR,
    rad.NEGATIVE_ANCHOR,
)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def charge_identities():
    """sigma enters Lambda, and the varied charge keeps the cubic factor."""
    c_theta, lhat, ell, number, zeta, r0, kappa, zed = s.symbols(
        'c_theta Lhat ell N zeta_s r0 kappa z', real=True,
    )
    pi = s.pi
    sigma = 2 * c_theta * lhat
    lambda_theta = (pi * ell**2 / number) * sigma
    lambda_zeta = pi * zeta * ell**2 / number
    coupled = lambda_zeta + lambda_theta
    closed = (pi * ell**2 / number) * (zeta + sigma)
    zero(s.together(coupled - closed))
    zero(s.together(lambda_theta - (2 * pi * ell**2 / number) * c_theta * lhat))
    factor = (r0 + kappa * zed) * (r0 + 3 * kappa * zed)
    zero(s.expand(factor - (r0**2 + 4 * kappa * r0 * zed + 3 * kappa**2 * zed**2)))
    return dict(
        sigma='2*c_theta*Lhat^{uu}',
        lambda_theta='(pi*ell**2/N)*sigma',
        Lambda='lambda_zeta+lambda_theta=(pi*ell**2/N)*(zeta_s+sigma)',
        charge_factor='(r0+kappa*z)*(r0+3*kappa*z)=r0**2+4*kappa*r0*z+3*kappa**2*z**2',
        c_theta_set_to_zero=False,
        sigma_set_to_zero=False,
    )


def bases(grid):
    stored = {}
    for z0, v0 in rad.corners():
        history = rad.march(rad.ETA, 0.0, 0.0, z0, v0, 'DOP853')
        if history is None:
            raise RuntimeError('base march failed')
        stored[(z0, v0)] = rad.sample_state(history, grid)[:2]
    return stored


def kappa_window(lam, grid, stored):
    positive = []
    negative = []
    for z0, v0 in rad.corners():
        base = stored[(z0, v0)]
        pos, pos_edge = rad.bracket_kappa(z0, v0, rad.ETA, lam, 1.0, grid, base)
        neg, neg_edge = rad.bracket_kappa(z0, v0, rad.ETA, lam, -1.0, grid, base)
        positive.append(dict(z0=z0, v0=v0, kappa=pos, edge=pos_edge))
        negative.append(dict(z0=z0, v0=v0, kappa=-neg, edge=neg_edge))
        print('kappa-corner', lam, z0, v0, -neg, pos, flush=True)
    return _shared_kappa(lam, positive, negative)


def _shared_kappa(lam, positive, negative):
    shared_positive = min(row['kappa'] for row in positive)
    shared_negative = max(row['kappa'] for row in negative)
    pos_worst = min(positive, key=lambda row: row['kappa'])
    neg_worst = max(negative, key=lambda row: row['kappa'])
    recorded = RECORDED_KAPPA.get(lam)
    if recorded is not None:
        assert abs(shared_negative - recorded[0]) < 1e-15, (lam, shared_negative, recorded[0])
        assert abs(shared_positive - recorded[1]) < 1e-15, (lam, shared_positive, recorded[1])
    print('kappa-sample', lam, shared_negative, shared_positive, flush=True)
    return dict(
        lambda_value=lam,
        kappa_negative=shared_negative,
        kappa_positive=shared_positive,
        positive_binding_budget=pos_worst['edge'].get('binding_budget'),
        negative_binding_budget=neg_worst['edge'].get('binding_budget'),
        positive_edge=_edge(pos_worst),
        negative_edge=_edge(neg_worst),
    )


def _edge(row):
    return dict(
        z0=row['z0'], v0=row['v0'],
        switch=row['edge'].get('sup_abs_delta_z_on_switch'),
        post_switch=row['edge'].get('post_switch_sup_abs_delta_z'),
        energy=row['edge'].get('abs_delta_e_f'),
        admissible=row['edge'].get('admissible'),
        binding_budget=row['edge'].get('binding_budget'),
    )


def bracket_lambda(z0, v0, kappa, sign, grid, base):
    lo = 0.0
    hi = 1.0
    while True:
        measured = rad.measure(
            rad.march(rad.ETA, sign * hi, kappa, z0, v0, 'DOP853'), base, grid, rad.ETA,
        )
        if not measured['admissible']:
            break
        hi *= 2.0
        if hi > 1.0e6:
            raise RuntimeError('joint Lambda budget remained admissible past 1e6')
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        measured = rad.measure(
            rad.march(rad.ETA, sign * mid, kappa, z0, v0, 'DOP853'), base, grid, rad.ETA,
        )
        if measured['admissible']:
            lo = mid
        else:
            hi = mid
    edge = rad.measure(
        rad.march(rad.ETA, sign * lo, kappa, z0, v0, 'DOP853'), base, grid, rad.ETA,
    )
    return lo, edge


def lambda_window(kappa, grid, stored):
    positive = []
    negative = []
    for z0, v0 in rad.corners():
        base = stored[(z0, v0)]
        pos, pos_edge = bracket_lambda(z0, v0, kappa, 1.0, grid, base)
        neg, neg_edge = bracket_lambda(z0, v0, kappa, -1.0, grid, base)
        positive.append(dict(z0=z0, v0=v0, value=pos, edge=pos_edge))
        negative.append(dict(z0=z0, v0=v0, value=-neg, edge=neg_edge))
        print('lambda-corner', kappa, z0, v0, -neg, pos, flush=True)
    shared_positive = min(row['value'] for row in positive)
    shared_negative = max(row['value'] for row in negative)
    pos_worst = min(positive, key=lambda row: row['value'])
    neg_worst = max(negative, key=lambda row: row['value'])
    print('lambda-sample', kappa, shared_negative, shared_positive, flush=True)
    return dict(
        kappa=kappa,
        lambda_negative=shared_negative,
        lambda_positive=shared_positive,
        positive_binding_budget=pos_worst['edge'].get('binding_budget'),
        negative_binding_budget=neg_worst['edge'].get('binding_budget'),
        positive_edge=_edge(pos_worst),
        negative_edge=_edge(neg_worst),
    )


def prescribed_reduction(grid, stored):
    rows = []
    for lam in (rad.POSITIVE_ANCHOR, rad.NEGATIVE_ANCHOR):
        corner_rows = []
        for z0, v0 in rad.corners():
            measured = rad.measure(
                rad.march(rad.ETA, lam, 0.0, z0, v0, 'DOP853'),
                stored[(z0, v0)], grid, rad.ETA,
            )
            assert measured['admissible'], (lam, z0, v0, measured)
            corner_rows.append(dict(
                z0=z0, v0=v0, admissible=True,
                binding_budget=measured['binding_budget'],
                switch=measured['sup_abs_delta_z_on_switch'],
                post_switch=measured['post_switch_sup_abs_delta_z'],
                energy=measured['abs_delta_e_f'],
            ))
            print('kappa0-anchor', lam, z0, v0, measured['binding_budget'], flush=True)
        rows.append(dict(lambda_value=lam, all_corners_admissible=True, corners=corner_rows))
    lam = rad.OVERSHOOT * rad.POSITIVE_ANCHOR
    overshoot = rad.measure(
        rad.march(rad.ETA, lam, 0.0, rad.WORST_Z0, rad.WORST_V0, 'DOP853'),
        stored[(rad.WORST_Z0, rad.WORST_V0)], grid, rad.ETA,
    )
    assert not overshoot['admissible'], overshoot
    print('kappa0-overshoot', lam, overshoot.get('binding_budget'), flush=True)
    return dict(
        anchors=rows,
        overshoot_lambda=lam,
        overshoot_corner=dict(z0=rad.WORST_Z0, v0=rad.WORST_V0),
        overshoot_kappa_zero_admissible=False,
        overshoot_binding_budget=overshoot.get('binding_budget'),
        overshoot_switch=overshoot.get('sup_abs_delta_z_on_switch'),
        overshoot_post_switch=overshoot.get('post_switch_sup_abs_delta_z'),
        overshoot_energy=overshoot.get('abs_delta_e_f'),
    )


def full_edge_probe(zero_sample, grid, stored):
    probes = []
    for edge, kappa in (
        (zero_sample['positive_edge'], zero_sample['kappa_positive']),
        (zero_sample['negative_edge'], zero_sample['kappa_negative']),
    ):
        z0 = edge['z0']
        v0 = edge['v0']
        samples = []
        for lam in PROBE_LAMBDAS:
            measured = rad.measure(
                rad.march(rad.ETA, lam, kappa, z0, v0, 'DOP853'),
                stored[(z0, v0)], grid, rad.ETA,
            )
            samples.append(dict(
                lambda_value=lam,
                admissible=measured['admissible'],
                binding_budget=measured.get('binding_budget'),
                switch=measured.get('sup_abs_delta_z_on_switch'),
                post_switch=measured.get('post_switch_sup_abs_delta_z'),
                energy=measured.get('abs_delta_e_f'),
            ))
            print('full-edge', kappa, lam, measured['admissible'], flush=True)
        probes.append(dict(kappa=kappa, corner=dict(z0=z0, v0=v0), samples=samples))
    return dict(
        probe_is_an_empty_interval_claim=False,
        edges=probes,
    )


def independent_edge(sample, grid):
    z0 = sample['positive_edge']['z0']
    v0 = sample['positive_edge']['v0']
    lam = sample['lambda_positive']
    kappa = sample['kappa']
    dop = rad.sample_state(rad.march(rad.ETA, lam, kappa, z0, v0, 'DOP853'), grid)[:2]
    other = rad.sample_state(rad.march(rad.ETA, lam, kappa, z0, v0, 'RK45'), grid)[:2]
    gap = float(np.max(np.abs(dop - other)))
    assert gap < 1e-7, gap
    print('independent', lam, kappa, gap, flush=True)
    return dict(
        lambda_value=lam, kappa=kappa, corner=dict(z0=z0, v0=v0), dop853_rk45_gap=gap,
    )


def main():
    warnings.filterwarnings('error', category=scipy.integrate.IntegrationWarning)
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    corner_list = rad.corners()
    assert len(corner_list) == 6, len(corner_list)
    primitive = alloc.shape_identities().pop('primitive_polynomial')
    force, _ = alloc.make_forces(alloc.antiderivative_factory(primitive))
    identities = charge_identities()
    drive = rad.drive_regression(force)
    anchor = rad.anchor_regression(force)
    print('regressions', drive['max_abs_drive_error'], anchor['anchor_trajectory_gap'], flush=True)
    grid = rad.grids()
    stored = bases(grid)
    reduction = prescribed_reduction(grid, stored)
    zero_kappa = kappa_window(0.0, grid, stored)
    half_kappa_ranges = [
        lambda_window(0.5 * zero_kappa['kappa_positive'], grid, stored),
        lambda_window(0.5 * zero_kappa['kappa_negative'], grid, stored),
    ]
    half_anchor_ranges = [
        kappa_window(0.5 * rad.POSITIVE_ANCHOR, grid, stored),
        kappa_window(0.5 * rad.NEGATIVE_ANCHOR, grid, stored),
    ]
    probe = full_edge_probe(zero_kappa, grid, stored)
    edge = independent_edge(half_kappa_ranges[0], grid)
    result = dict(
        status='PASS',
        verdict='finite_weyl_product_on_realized_angle',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        area_budget='not evaluated; no replacement target was written before the test',
        c_theta_set_to_zero=False,
        sigma_set_to_zero=False,
        corner_count=len(corner_list),
        identified_input=(
            'sigma=2*c_theta*Lhat^{uu} enters as lambda_theta inside '
            'Lambda=lambda_zeta+lambda_theta on the realized delay; '
            'zeta_s stays unmatched inside that sum; the O(epsilon) '
            'product-throat weight is a geometric factor'
        ),
        identities=identities,
        drive_regression=drive,
        anchor_regression=anchor,
        prescribed_reduction=reduction,
        lambda_zero_kappa=zero_kappa,
        lambda_ranges_at_half_kappa=half_kappa_ranges,
        kappa_ranges_at_half_anchors=half_anchor_ranges,
        full_edge_probe=probe,
        independent_integrator=edge,
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
        corner_count=result['corner_count'],
        lambda_at_half_positive_kappa=[
            half_kappa_ranges[0]['lambda_negative'],
            half_kappa_ranges[0]['lambda_positive'],
        ],
        lambda_at_half_negative_kappa=[
            half_kappa_ranges[1]['lambda_negative'],
            half_kappa_ranges[1]['lambda_positive'],
        ],
        kappa_at_half_positive_anchor=[
            half_anchor_ranges[0]['kappa_negative'],
            half_anchor_ranges[0]['kappa_positive'],
        ],
        kappa_at_half_negative_anchor=[
            half_anchor_ranges[1]['kappa_negative'],
            half_anchor_ranges[1]['kappa_positive'],
        ],
        independent_gap=edge['dop853_rk45_gap'],
        s074_gate_closed=False,
        opening_achieved=False,
        c_theta_set_to_zero=False,
    )
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
