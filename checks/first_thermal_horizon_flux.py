"""Audit the first thermal full-tower horizon constraint, not stability."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.special import expit
from numpy.polynomial.legendre import leggauss
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = 'first-thermal-horizon-flux'
OWN = ['checks/first_thermal_horizon_flux.py', 'notes/'+STEM+'.md',
       'evidence/first-thermal-horizon-source-review.json']
PREVIOUS = 'causal-horizon-and-zero-mode-response'
ROWS, NUMERICS = [], []
PROTOCOL = {
    'scope': 'One exact full Dirac horizon flux value for separated compact input and fixed HH preparation',
    'matrix_fixtures': [[1.7, 1, False], [3.1, 1, False],
                        [1.7, 2, False], [2.3, 1, True]],
    'quadrature_nodes': [48, 96], 'matrix_absolute_error_limit': 2e-11,
    'tail_fixtures': [[4., 1., 1., 0, 64], [100., 10., .8, 4, 64],
                      [float(2000**4), float(2000**4)/350., 1., 8, 16]],
    'tail_cutoffs': [1024, 4096],
    'full_quantum_kernel_evaluated': False, 'formal_verification': False,
}


def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    row = {'name': name, 'passed': bool(ok), 'detail': str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    value = S.simplify(S.sympify(expression).doit())
    check(name, value == 0, value)


def prior():
    return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]


def paths():
    return sorted(set(OWN) | set(prior()['input_sha256_after']) |
                  {'evidence/'+PREVIOUS+'.json'})


def provenance():
    p = prior()
    check('predecessor_passed', p['status'] == 'passed' and len(p['checks']) == 98
          and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable', p['input_sha256_before'] == p['input_sha256_after'])
    bad = [f for f, h in p['input_sha256_after'].items()
           if sha256((ROOT/f).read_bytes()).hexdigest() != h]
    check('predecessor_hashes', not bad, bad)
    review = json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned', all(f in p['input_sha256_after']
                                      for f in review['retained_inputs']))
    check('source_scope', len(review['primary_sources']) == 3
          and review['full_angular_tower_retained']
          and review['no_optical_mass_gap_used']
          and review['compact_source_and_fixed_state_required']
          and not review['pole_skipping_claimed'])


def exact_checks():
    U, V, x, k, xr = S.symbols('U V x k xr', nonzero=True)
    A, B, C = S.symbols('A B C')
    coord = [U, V]
    xi = S.Matrix([-k*U, k*V])
    a, b, c = (S.Function(s)(-U*V) for s in ['a', 'b', 'c'])
    tensor = S.Matrix([[V**3*a, V*b], [V*b, U*c]])
    for i, j in [(0, 0), (0, 1), (1, 1)]:
        lie = sum(xi[l]*S.diff(tensor[i, j], coord[l])
                  + tensor[l, j]*S.diff(xi[l], coord[i])
                  + tensor[i, l]*S.diff(xi[l], coord[j]) for l in range(2))
        zero('first_harmonic_Lie_'+str(i)+str(j), lie-k*tensor[i, j])
    J = S.Matrix([[-k*U, -xr/V], [k*V, 0]])
    T = S.Matrix([[V**3*A, V*B], [V*B, U*C]])
    EF = S.simplify((J.T*T*J).subs(U, -x/V)/V)
    expected = S.Matrix([[k*k*x*(x*A+2*B-C), -k*xr*(x*A+B)],
                         [-k*xr*(x*A+B), xr*xr*A]])
    for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        zero('covariant_transform_'+str(i)+str(j), EF[i, j]-expected[i, j])
    metric = S.Matrix([[-2*k*x/xr, 1], [1, 0]])
    mixed = metric.inv()*EF
    zero('mixed_flux_raised_before_limit', mixed[1, 0]+k*k*x*(C+x*A))
    zero('first_harmonic_horizon_flux', S.limit(mixed[1, 0], x, 0))
    zero('horizon_pressure_retained', S.limit(EF[1, 1], x, 0)-xr*xr*A)
    zero('horizon_density_retained', S.limit(mixed[0, 0], x, 0)+k*xr*B)
    # Evaluate the off-source spherical energy Ward identity independently.
    r, rh, aa = S.symbols('r rh aa', positive=True)
    xx = aa*(r-rh)
    flux = -k*k*xx*(C+xx*A)
    density = -k*aa*(xx*A+B)
    ward = (S.diff(r*r*flux, r)+k*r*r*density).subs(r, rh)
    zero('horizon_Ward_coefficient', ward+aa*k*k*rh*rh*(B+C))
    zero('horizon_Ward_relation', ward.subs(C, -B))
    # Independent full divergence in the flat two-dimensional product control.
    flat = S.Matrix([[0, -1], [-1, 0]])
    nonzero = S.Matrix([[0, V], [V, -U]])
    raised = flat.inv()*nonzero
    for j in range(2):
        zero('nonzero_smooth_conserved_control_'+str(j),
             sum(S.diff(raised[i, j], coord[i]) for i in range(2)))
    zero('nonzero_trace_control', S.trace(raised)+2*V)
    # Future-regular null radiation is singular at the bifurcation surface.
    amplitude = V/(4*S.pi*r*r)
    future_VV = amplitude/(k*k*V*V)
    zero('future_regular_bifurcation_singularity', V*future_VV-1/(4*S.pi*k*k*r*r))
    zero('future_regular_null_energy_Ward', S.diff(r*r*amplitude, r))
    check('future_regular_flux_not_zero', amplitude.subs({V: 1, r: 1}) != 0)
    # Smooth n=2 tensors permit nonzero horizon flux.
    T2 = S.Matrix([[V**4*A, V**2*B], [V**2*B, C]])
    EF2 = S.simplify((J.T*T2*J).subs(U, -x/V)/(V*V))
    flux2 = (metric.inv()*EF2)[1, 0]
    zero('second_harmonic_mixed_flux', flux2-k*k*(C-x*x*A))
    zero('second_harmonic_nonzero_horizon', S.limit(flux2, x, 0)-k*k*C)
    # Explicit angular integral and monotonicity sufficient condition.
    q, bb, z = S.symbols('q bb z', positive=True)
    tail = 4*S.exp(-bb*z)*(q*(z/bb+1/bb**2)
                            + 2*(z**3/bb+3*z*z/bb**2+6*z/bb**3+6/bb**4))
    zero('tail_antiderivative', S.diff(tail, z)+4*(q*z+2*z**3)*S.exp(-bb*z))
    zero('tail_at_infinity', S.limit(tail, z, S.oo))
    zero('decreasing_envelope_margin', 1/z-2/(q+2*z)-q/(z*(q+2*z)))
    n, CC = S.symbols('n CC', positive=True)
    zero('proper_mass_lower_bound_remainder', n*(n+q)/CC-q*n/CC-n*n/CC)
    # The exact LLL factor is consistent with the full-tower constraint.
    g = S.symbols('g')
    zero('previous_LLL_frequency_zero', (g*(k*k-g*g)).subs(g, k))


def matrix_fixtures():
    H0 = np.array([[-.7, .2+.1j, 0], [.2-.1j, .1, .3j], [0, -.3j, .8]])
    O = np.array([[.2, .1j, .3], [-.1j, -.4, .2], [.3, .2, .7]])
    B = np.array([[-.1, .4, .1j], [.4, .3, .2], [-.1j, .2, .2]])
    W = np.diag([.3, -.2, .4])
    for index, (bt, harmonic, degenerate) in enumerate(PROTOCOL['matrix_fixtures']):
        H = np.diag([-.4, -.4, .7]) if degenerate else H0
        e, P = np.linalg.eigh(H)
        cov = (P*expit(-bt*e))@P.conj().T
        gamma = 2*np.pi*harmonic/bt
        generator = gamma*np.eye(9)+1j*(np.kron(np.eye(3), H)-np.kron(H.T, np.eye(3)))
        forcing = -1j*(B@cov-cov@B)
        delta = np.linalg.solve(generator, forcing.reshape(9, order='F')).reshape((3, 3), order='F')
        contact = np.trace(cov@W)
        causal = np.trace(O@delta)+contact
        euclidean = []
        for nodes in PROTOCOL['quadrature_nodes']:
            zz, ww = leggauss(nodes)
            integral = 0j
            for tau, weight in zip((zz+1)*bt/2, ww*bt/2):
                kernel = np.trace(cov@expm(tau*H)@O@expm(-tau*H)@(np.eye(3)-cov)@B)
                integral += weight*np.exp(1j*gamma*tau)*kernel
            euclidean.append(contact-integral)
        error = abs(causal-euclidean[-1])
        refinement = abs(euclidean[-1]-euclidean[0])
        limit = PROTOCOL['matrix_absolute_error_limit']
        check('KMS_resolvent_vs_Euclidean_'+str(index), error < limit, error)
        check('KMS_quadrature_refinement_'+str(index), refinement < limit, refinement)
        check('KMS_dropped_contact_detected_'+str(index),
              abs(causal-(euclidean[-1]-contact)) > .05, contact)
        check('KMS_noncommuting_probe_'+str(index), np.linalg.norm(H@B-B@H) > .1)
        NUMERICS.append({'fixture': 'KMS_matrix', 'beta': bt, 'harmonic': harmonic,
                         'degenerate': degenerate, 'causal': [causal.real, causal.imag],
                         'euclidean': [[v.real, v.imag] for v in euclidean],
                         'contact': [contact.real, contact.imag],
                         'absolute_error': error, 'refinement': refinement})


def tail_bound(q, b, N):
    z = np.sqrt(N)
    return 4*np.exp(-b*z)*(q*(z/b+1/b**2)
                           + 2*(z**3/b+3*N/b**2+6*z/b**3+6/b**4))


def tail_fixtures():
    for index, (q, C, t, p, N) in enumerate(PROTOCOL['tail_fixtures']):
        b = t*np.sqrt(q/C)/2
        K = 1. if p == 0 else (1+2*p/t)**p
        upper = K*tail_bound(q, b, N)
        check('tail_monotone_range_'+str(index), N >= max(1., 4/(b*b)))
        values = []
        for cutoff in PROTOCOL['tail_cutoffs']:
            n = np.arange(N+1, cutoff+1, dtype=float)
            M = np.sqrt(n*(n+q)/C)
            actual = 2*(q+2*n)*np.exp(p*np.log1p(M)-t*M)
            envelope = K*2*(q+2*n)*np.exp(-b*np.sqrt(n))
            check('pointwise_tail_majorant_'+str(index)+'_'+str(cutoff),
                  bool(np.all(actual <= envelope*(1+1e-13))))
            remainder = K*tail_bound(q, b, cutoff)
            envelope_with_remainder = float(np.sum(envelope))+remainder
            actual_with_remainder = float(np.sum(actual))+remainder
            check('integral_tail_majorant_'+str(index)+'_'+str(cutoff),
                  envelope_with_remainder <= upper*(1+1e-13),
                  envelope_with_remainder/upper)
            values.append({'cutoff': cutoff, 'actual_partial': float(np.sum(actual)),
                           'envelope_partial': float(np.sum(envelope)),
                           'analytic_remainder': remainder,
                           'actual_plus_remainder': actual_with_remainder,
                           'envelope_plus_remainder': envelope_with_remainder})
        check('tail_remainder_refines_'+str(index),
              values[1]['analytic_remainder'] <= values[0]['analytic_remainder'])
        NUMERICS.append({'fixture': 'angular_tail_envelope', 'q': q, 'C': C,
                         't': t, 'p': p, 'N': N, 'b': b, 'K': K,
                         'analytic_upper': upper, 'cutoffs': values,
                         'physical_kernel_error_bar': False})


def main():
    sys.stdout.reconfigure(encoding='utf-8', newline='\n')
    sys.stderr.reconfigure(encoding='utf-8', newline='\n')
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before, error = hashes(), None
    print(json.dumps({'event': 'start', 'story': 'S168 first thermal horizon flux',
                      'input_count': len(ps)}), flush=True)
    try:
        provenance()
        exact_checks()
        matrix_fixtures()
        tail_fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r['passed'] for r in ROWS)
    run = {'story': 'S168 first thermal horizon flux',
           'created_utc': datetime.now(timezone.utc).isoformat(),
           'status': 'passed' if good else 'failed', 'exception': error,
           'python': sys.version, 'sympy': S.__version__, 'numpy': np.__version__,
           'scipy': scipy.__version__, 'protocol': PROTOCOL,
           'input_sha256_before': before, 'input_sha256_after': after,
           'unchanged_inputs': before == after, 'checks': ROWS, 'numerics': NUMERICS,
           'first_thermal_full_Dirac_horizon_flux_zero_derived': good,
           'separated_massive_sum_smoothness_derived': good,
           'fixed_contacts_retained': good,
           'full_quantum_kernel_evaluated': False, 'general_full_horizon_trace_proved': False,
           'uniform_infinite_time_bound_proved': False, 'physical_mass_return_proved': False,
           'massive_tower_discarded': False, 'optical_mass_gap_assumed': False,
           'counterterm_changed': False, 'physical_metric_boundary_domain_selected': False,
           'global_shape_inverse_proved': False, 'actual_metric_mode_constructed': False,
           'actual_charged_stability_proved': False, 'actual_wormhole_instability_proved': False,
           'pole_skipping_claimed': False, 'retained_fast_factor_removed': False,
           'order_reduction_adopted': False, 'S167_complete': False,
           'S168_complete': False, 'S077_complete': False, 'goal_complete': False}
    if not good:
        for kind, path in [('checker', OWN[0]), ('note', OWN[1]), ('review', OWN[2])]:
            run['failed_'+kind+'_source_utf8'] = (ROOT/path).read_text(encoding='utf-8')
    dest = ROOT/'evidence'/(STEM+'.json')
    data = json.loads(dest.read_bytes()) if dest.exists() else {'runs': []}
    data['runs'].append(run)
    dest.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': run['status'], 'checks': len(ROWS),
                      'unchanged_inputs': before == after}), flush=True)
    return 0 if good else 1


if __name__ == '__main__':
    raise SystemExit(main())
