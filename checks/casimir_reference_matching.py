"""Independent first-order Einstein-Maxwell contraction for a magnetic throat.

Reproduces the invariant-vacuum plus leading LLL Casimir baseline. The
supplied state difference is an assumption, not a full global determinant.
No Lean invocation and no higher-order or physical stability certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def geometry(metric, coordinates):
    inverse = metric.inv()
    gamma = {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                gamma[a, b, c] = sp.simplify(sum(
                    inverse[a, d] * (sp.diff(metric[d, c], coordinates[b])
                    + sp.diff(metric[d, b], coordinates[c])
                    - sp.diff(metric[b, c], coordinates[d])) / 2 for d in range(4)))
    ricci = sp.Matrix(4, 4, lambda b, d: sp.simplify(sum(
        sp.diff(gamma[a, b, d], coordinates[a])
        - sp.diff(gamma[a, a, b], coordinates[d])
        + sum(gamma[a, a, e] * gamma[e, b, d]
              - gamma[a, d, e] * gamma[e, a, b] for e in range(4))
        for a in range(4))))
    scalar = sp.simplify(sp.trace(inverse * ricci))
    return inverse, gamma, ricci - metric * scalar / 2, scalar


def verify():
    checks = {}

    def zero(name, value):
        reduced = sp.simplify(value)
        assert reduced == 0, (name, reduced)
        checks[name] = 'exact_zero'

    t, r, u, phi = sp.symbols('tau rho u phi', real=True)
    eps, U, V, chi = sp.symbols('epsilon U V chi', real=True)
    q, re, G = sp.symbols('q r_e G', positive=True)
    coords = [t, r, u, phi]
    f, sphere = 1 + r**2, 1 - u**2
    S, C = sp.Function('S')(r), sp.Function('C')(r)
    metric = sp.diag(-S, 1/S, C/sphere, C*sphere)
    inverse, gamma, einstein, scalar = geometry(metric, coords)
    # Work in units r_e=1; kappa/g_c^2=8/q^2 follows from fixed magnetic charge.
    field = sp.zeros(4)
    field[2, 3], field[3, 2] = q/2, -q/2
    upper_field = inverse * field * inverse
    field_squared = sum(field[i, j]*upper_field[i, j]
                        for i in range(4) for j in range(4))
    magnetic = (8/q**2) * (field * inverse * field.T - metric * field_squared/4)
    residual = einstein - magnetic
    normalizers = [1/S, S, sphere/C, 1/(sphere*C)]
    c1, s1 = sp.Function('c')(r), sp.Function('s')(r)
    substitution = {S: f + eps*s1, C: 1 + eps*c1}

    def coefficient(expression):
        expanded = expression.subs(substitution, simultaneous=True).doit()
        return sp.simplify(sp.diff(expanded, eps).subs(eps, 0))

    predicted = [c1 - f*sp.diff(c1, r, 2) - r*sp.diff(c1, r),
                 -c1 + r*sp.diff(c1, r),
                 sp.diff(s1, r, 2)/2 + f*sp.diff(c1, r, 2)/2
                 + r*sp.diff(c1, r) + 2*c1]
    predicted.append(predicted[2])
    linear_components = []
    for i in range(4):
        expression = sp.simplify(normalizers[i] * residual[i, i])
        zero(f'classical_Einstein_Maxwell_{i}',
             expression.subs({S: f, C: 1}, simultaneous=True).doit())
        linear_components.append(coefficient(expression))
        zero(f'independent_linearized_component_{i}', linear_components[i]-predicted[i])
    for i in range(4):
        for j in range(i+1, 4):
            zero(f'exact_off_diagonal_{i}_{j}', residual[i, j])
    # sqrt(-det(g))=C on the positive-radius coordinate patch.
    for nu in range(4):
        zero(f'exact_magnetic_Maxwell_{nu}',
             sum(sp.diff(C*upper_field[mu, nu], coords[mu]) for mu in range(4))/C)

    c_solution = chi*(1 + r*sp.atan(r)) - U
    casimir_s = -r**2 - r*(r**2+3)*sp.atan(r) + sp.log(f)
    s_solution = (V+2*U)*r**2 + chi*casimir_s
    solution = {c1: c_solution, s1: s_solution}
    source_normalized = [-U-chi/f, U-chi/f, V, V]
    for i in range(4):
        zero(f'combined_source_equation_{i}',
             linear_components[i].subs(solution, simultaneous=True).doit()-source_normalized[i])
    zero('area_opening_equation', sp.diff(c_solution, r, 2)-2*chi/f**2)
    zero('Casimir_lapse_second_derivative', sp.diff(casimir_s, r, 2)+6*(1+r*sp.atan(r)))
    zero('reflection_c_prime', sp.diff(c_solution, r).subs(r, 0))
    zero('reflection_s_prime', sp.diff(s_solution, r).subs(r, 0))
    zero('central_time_normalization', s_solution.subs(r, 0))

    background = {S: f, C: 1}
    inverse0 = inverse.subs(background, simultaneous=True)
    gamma0 = {key: value.subs(background, simultaneous=True).doit()
              for key, value in gamma.items()}
    source = sp.diag(-U*f-chi, U/f-chi/f**2, V/sphere, V*sphere)
    mixed_source = inverse0 * source
    for nu in range(4):
        divergence = sum(sp.diff(mixed_source[mu, nu], coords[mu])
                         + sum(gamma0[mu, mu, lam]*mixed_source[lam, nu]
                               - gamma0[lam, mu, nu]*mixed_source[mu, lam]
                               for lam in range(4)) for mu in range(4))
        zero(f'background_source_conservation_{nu}', divergence)
    zero('source_trace', sp.trace(mixed_source)-2*(U+V))
    zero('curvature_trace_equation',
         coefficient(scalar).subs(solution, simultaneous=True).doit()+2*(U+V))
    null_vector = sp.Matrix([1/f, 1, 0, 0])
    zero('background_radial_null', (null_vector.T*metric.subs(background)*null_vector)[0])
    null_source = (null_vector.T*source*null_vector)[0]
    zero('radial_null_source', null_source+2*chi/f**2)
    zero('vacuum_null_source', null_source.subs(chi, 0))
    first_radius = re*c_solution/2
    zero('areal_radius_second_derivative_coefficient', sp.diff(first_radius, r, 2)-re*chi/f**2)
    zero('overlap_radius_slope_coefficient', sp.limit(sp.diff(first_radius, r), r, sp.oo)-sp.pi*re*chi/4)

    # The vacuum-only metric is exactly a product after rescaling both AdS coordinates.
    A, rh = sp.symbols('A rho_hat', positive=True)
    zero('reference_AdS_time_rescaling', (-(1+A*r*r)/A).subs(r, rh/sp.sqrt(A))+(1+rh**2)/A)
    zero('reference_AdS_radial_rescaling', (1/(A*(1+A*r*r))).subs(r, rh/sp.sqrt(A))-1/(A*(1+rh**2)))
    radius_ratio = 1/((1+eps*(V+2*U))*(1-eps*U))
    zero('reference_radius_ratio_shift', sp.diff(radius_ratio, eps).subs(eps, 0)+U+V)
    zero('W22_area_specialization', c_solution.subs({U: -sp.Rational(1, 3), V: sp.Rational(1, 3), chi: 1})
         -r*sp.atan(r)-sp.Rational(4, 3))
    zero('W22_lapse_specialization', s_solution.subs({U: -sp.Rational(1, 3), V: sp.Rational(1, 3), chi: 1})
         +sp.Rational(4, 3)*r**2+r*(r**2+3)*sp.atan(r)-sp.log(f))
    alpha = G*q/(4*sp.pi*re**2)
    leading_length = 4*re/(sp.pi*alpha)
    zero('MMP_leading_length', leading_length-16*re**3/(G*q))
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = verify()
    result = {
        'status': 'passed_exact_first_order_baseline_checks',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'checks': checks,
        'identities_checked': len(checks),
        'metric_contracted_independently': True,
        'source': 'Specified invariant pressures plus leading long-throat LLL state difference',
        'order': 'First order in a common loop bookkeeping parameter at fixed charge on a finite patch',
        'known_baseline_reproduced': True,
        'global_quantum_state_matched': False,
        'higher_order_response_controlled': False,
        'corrected_global_transit_time_computed': False,
        'physical_stability_proved': False,
        'novel_physics_result': False,
        'lean_certified': False,
        'lean_invocations': 0
    }
    args.output.write_bytes((json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps({'status': result['status'], 'identities_checked': len(checks),
                      'script_sha256': result['script_sha256']}))


if __name__ == '__main__':
    main()
