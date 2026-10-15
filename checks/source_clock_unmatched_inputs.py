"""Classify the four unmatched inputs on the retained rule-D clock.

The comparison for each input is in notes/source-clock-unmatched-inputs.md.
This checker does not re-bracket Lambda and does not edit the joint,
radiation, controller, or allocation checkers.
"""
import argparse
import hashlib
import inspect
import json
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import sympy as s

import source_clock_allocation as alloc
import source_clock_joint as joint
import source_clock_radiation as rad


ROOT = Path(__file__).resolve().parents[1]
JOINT_RECEIPT = ROOT / 'evidence' / 'source-clock-joint.json'


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def recorded_point():
    """Use one admissible S084 point. Do not compute a new interval."""
    receipt = json.loads(JOINT_RECEIPT.read_text(encoding='utf-8'))
    sample = receipt['lambda_ranges_at_half_kappa'][0]
    edge = sample['positive_edge']
    return dict(
        kappa=sample['kappa'],
        lambda_positive=sample['lambda_positive'],
        lambda_negative=sample['lambda_negative'],
        z0=edge['z0'],
        v0=edge['v0'],
        recorded_lambda_intervals=[
            [row['lambda_negative'], row['lambda_positive']]
            for row in receipt['lambda_ranges_at_half_kappa']
        ],
    )


def shift_identity():
    """Opposite shifts leave Lambda fixed. Neither coefficient is zero."""
    declared = joint.charge_identities()
    assert declared['Lambda'] == (
        'lambda_zeta+lambda_theta=(pi*ell**2/N)*(zeta_s+sigma)'
    )
    assert declared['c_theta_set_to_zero'] is False
    assert declared['sigma_set_to_zero'] is False
    zeta, sigma, shift, ell, number = s.symbols(
        'zeta_s sigma a ell N', real=True,
    )
    pi = s.pi
    coupled = (pi * ell**2 / number) * (zeta + sigma)
    shifted = (pi * ell**2 / number) * ((zeta + shift) + (sigma - shift))
    zero(s.together(shifted - coupled))
    assert s.simplify(s.diff(coupled, zeta) - s.diff(coupled, sigma)) == 0
    assert s.simplify(s.diff(coupled, zeta)) != 0
    lam_zeta, lam_theta = s.symbols('lambda_zeta lambda_theta', real=True)
    total = lam_zeta + lam_theta
    zero(s.expand((lam_zeta + shift) + (lam_theta - shift) - total))
    return dict(
        identity='Lambda=(pi*ell**2/N)*(zeta_s+sigma)',
        shift='zeta_s+a with sigma-a leaves Lambda fixed',
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
    )


def clock_polynomial(lam, kappa, r0, zed, drive, extra):
    """Local acceleration. extra collects symbols the clock does not use."""
    return (
        -(1 + lam * r0**2) * zed
        - 4 * lam * kappa * r0 * zed**2
        - 3 * lam * kappa**2 * zed**3
        + s.Rational(1, 2) * drive
        + 0 * extra
    )


def absent_derivatives():
    lam, kappa, r0, zed, drive = s.symbols(
        'Lambda kappa r0 z drive', real=True,
    )
    collar, logdet, cj = s.symbols('collar_difference logdet cJ', real=True)
    extra = collar + logdet + cj
    accel = clock_polynomial(lam, kappa, r0, zed, drive, extra)
    # The added zero term must not change the polynomial that march uses.
    zero(s.expand(accel - clock_polynomial(lam, kappa, r0, zed, drive, 0)))
    for symbol in (collar, logdet, cj):
        zero(s.diff(accel, symbol))
    moved = s.diff(accel, lam)
    assert s.simplify(moved) != 0
    # cJ's operator is not the clock contact.
    jay, rate, kem = s.symbols('J rate kappa_em', real=True)
    delta = cj * kem**4 * jay**2 / 16
    zero(delta.subs(jay, 0))
    zero(s.diff(delta, jay).subs(jay, 0))
    assert cj in delta.free_symbols
    contact = rate**2
    assert jay in s.expand(delta).free_symbols
    assert rate not in s.expand(delta).free_symbols
    assert jay not in s.expand(contact).free_symbols
    return dict(
        collar_derivative=0,
        determinant_derivative=0,
        cJ_derivative=0,
        lambda_derivative_is_zero=False,
        magnetic_background_measures_cJ=False,
        cJ_set_to_zero=False,
        conditional_axial_window_adopted=False,
    )


def _march_or_fail(lam, kappa, z0, v0, method):
    history = rad.march(rad.ETA, lam, kappa, z0, v0, method)
    if history is None:
        raise RuntimeError('shipped march failed')
    return history


def state_gap(left, right, grid):
    a = rad.sample_state(left, grid)
    b = rad.sample_state(right, grid)
    return float(np.max(np.abs(a - b)))


def reconstructed_gap(history, lam, kappa):
    """Finite-difference acceleration against the local polynomial."""
    gaps = []
    eta = rad.ETA
    step = 1e-5
    for t0, t1, sol in history:
        if t1 - t0 <= 4 * step:
            continue
        times = np.linspace(t0 + 2 * step, t1 - 2 * step, 7)
        for time in times:
            now = sol(float(time))
            later = sol(float(time) + step)
            numeric = (later[1] - now[1]) / step
            zed = float(now[0])
            theta = float(now[2])
            charge = float(now[3])
            r0 = 4.0 * eta * float(alloc.angle_rate_factor(float(time)))
            drive = (theta - charge) / (4.0 * np.pi * eta)
            formula = (
                -(1.0 + lam * r0 * r0) * zed
                - 4.0 * lam * kappa * r0 * zed * zed
                - 3.0 * lam * kappa * kappa * zed * zed * zed
                + 0.5 * drive
            )
            gaps.append(abs(numeric - formula))
    if not gaps:
        raise RuntimeError('no reconstruction samples')
    return float(max(gaps))


def numerical(point):
    signature = inspect.signature(rad.march)
    assert list(signature.parameters) == ['eta', 'lam', 'kappa', 'z0', 'v0', 'method']
    lam = float(point['lambda_positive'])
    other = float(point['lambda_negative'])
    kappa = float(point['kappa'])
    z0 = float(point['z0'])
    v0 = float(point['v0'])
    # Both pieces nonzero. The sums are the recorded Lambda, not a new fit.
    assert 0.25 * lam != 0.0 and 0.75 * lam != 0.0
    same_left = _march_or_fail(0.25 * lam + 0.75 * lam, kappa, z0, v0, 'DOP853')
    same_right = _march_or_fail(0.75 * lam + 0.25 * lam, kappa, z0, v0, 'DOP853')
    different = _march_or_fail(other, kappa, z0, v0, 'DOP853')
    grid = np.linspace(0.0, alloc.T_END, 241)
    grid[-1] = alloc.T_END
    same_gap = state_gap(same_left, same_right, grid)
    sum_gap = state_gap(same_left, different, grid)
    assert same_gap < 1e-12, same_gap
    assert sum_gap > same_gap
    poly_gap = reconstructed_gap(same_left, lam, kappa)
    assert poly_gap < 1e-4, poly_gap
    return dict(
        kappa=kappa,
        lambda_positive=lam,
        lambda_negative=other,
        z0=z0,
        v0=v0,
        split_state_gap=same_gap,
        distinct_sum_state_gap=sum_gap,
        reconstruction_gap=poly_gap,
        march_parameters=list(signature.parameters),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    point = recorded_point()
    direct = shift_identity()
    absent = absent_derivatives()
    measured = numerical(point)
    result = dict(
        status='PASS',
        verdict='scoped_obstructions_on_retained_rule_D_clock',
        opening_achieved=False,
        s074_gate_closed=False,
        not_an_impossibility_theorem=True,
        area_budget_evaluated=False,
        area_budget='not evaluated; no replacement target was written before the test',
        zeta_s_set_to_zero=False,
        sigma_set_to_zero=False,
        cJ_set_to_zero=False,
        new_lambda_interval=False,
        recorded_lambda_intervals=point['recorded_lambda_intervals'],
        tasks=dict(
            zeta_s=dict(
                classification='unavailable physical input',
                admissible_range=None,
                identified_physical_input=False,
                split_state_gap=measured['split_state_gap'],
                distinct_sum_state_gap=measured['distinct_sum_state_gap'],
            ),
            unequal_collars=dict(
                classification='unresolved proof',
                admissible_range=None,
                identified_physical_input=False,
            ),
            nonlocal_interacting_determinant=dict(
                classification='unresolved proof',
                admissible_range=None,
                identified_physical_input=False,
            ),
            cJ=dict(
                classification='unavailable physical input',
                admissible_range=None,
                identified_physical_input=False,
                conditional_axial_window_adopted=False,
            ),
        ),
        direct_variation=direct,
        absent_derivatives=absent,
        numerical=measured,
        versions=dict(
            python=platform.python_version(),
            sympy=s.__version__,
            numpy=np.__version__,
        ),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(
        status=result['status'],
        split_state_gap=measured['split_state_gap'],
        distinct_sum_state_gap=measured['distinct_sum_state_gap'],
        reconstruction_gap=measured['reconstruction_gap'],
        opening_achieved=False,
        s074_gate_closed=False,
    )), flush=True)


if __name__ == '__main__':
    main()
