"""S132: ordered det_3 response; fixed-domain source model only, no Lean."""
import os
for _thread_name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_thread_name] = '1'
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import block_diag
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nonlinear_current_determinant as prior

ROOT = Path(__file__).resolve().parents[1]
TOL = dict(rounding=5e-11, derivative=2e-6, hessian=2e-5,
           inequality_relative=5e-10, noncommuting_minimum=1e-6)
STEPS = [2e-4, 1e-4]
DEGREES = [3, 6, 12]
RNG = np.random.default_rng(20261015)


def sha(name):
    return hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def norms(a):
    sv = np.linalg.svd(a, compute_uv=False)
    return float(np.sum(sv**3)**(1/3)), float(sv[0])


def value(k):
    # Independent finite determinant, including the two exact subtractions.
    # Unwrap the determinant phase using eigenvalues of I+K in the right
    # half-plane; this gives the logarithm reached from K=0 in ||K||op<1.
    sign, logabs = np.linalg.slogdet(np.eye(len(k))+k)
    principal = float(np.angle(sign))
    continued_phase = float(np.sum(np.angle(1+np.linalg.eigvals(k))))
    phase = principal+2*np.pi*round((continued_phase-principal)/(2*np.pi))
    return logabs+1j*phase-np.trace(k)+np.trace(k@k)/2


def first(k, h):
    r = np.linalg.solve(np.eye(len(k))+k, np.eye(len(k)))
    return np.trace(k@k@r@h)


def second(k, j, h):
    r = np.linalg.solve(np.eye(len(k))+k, np.eye(len(k)))
    return np.trace(((k@j+j@k)@r-k@k@r@j@r)@h)


def series(k, h, j, degree):
    p = [np.eye(len(k), dtype=complex)]
    for _ in range(degree):
        p.append(p[-1]@k)
    f, d, d2 = 0j, 0j, 0j
    for n in range(3, degree+1):
        sign = (-1)**(n+1)
        f += sign*np.trace(p[n])/n
        d += sign*np.trace(p[n-1]@h)
        d2 += sign*sum(np.trace(p[i]@j@p[n-2-i]@h) for i in range(n-1))
    return f, d, d2


def assert_bound(actual, ceiling):
    assert actual <= ceiling*(1+TOL['inequality_relative'])+TOL['rounding'], (actual, ceiling)


def matrix_case(k, h, j, label):
    s, a = norms(k)
    hh, _ = norms(h)
    jj, _ = norms(j)
    assert a < 1
    f, d, d2 = value(k), first(k, h), second(k, j, h)
    bound_d = s*s*hh/(1-a)
    bound_d2 = s*(2-a)*hh*jj/(1-a)**2
    assert_bound(abs(d), bound_d)
    assert_bound(abs(d2), bound_d2)
    symmetry_error = abs(d2-second(k, h, j))
    assert symmetry_error < TOL['rounding']
    commutator = float(np.linalg.norm(k@h-h@k))
    assert commutator > TOL['noncommuting_minimum']
    differences = []
    for step in STEPS:
        df = (value(k+step*h)-value(k-step*h))/(2*step)
        mixed = (value(k+step*h+step*j)-value(k+step*h-step*j)
                 -value(k-step*h+step*j)+value(k-step*h-step*j))/(4*step**2)
        # The path has K''=J, testing the extra DF[K''] term in (4).
        curved = (value(k+step*h+step**2*j/2)-2*f
                  +value(k-step*h+step**2*j/2))/step**2
        path_exact = second(k, h, h)+first(k, j)
        errors = dict(step=step, first=float(abs(df-d)),
                      mixed=float(abs(mixed-d2)), path_second=float(abs(curved-path_exact)))
        assert errors['first'] < TOL['derivative'], errors
        assert max(errors['mixed'], errors['path_second']) < TOL['hessian'], errors
        differences.append(errors)
    tails = []
    for degree in DEGREES:
        f_n, d_n, d2_n = series(k, h, j, degree)
        common = a**(degree-2)
        bounds = [s**3*common/((degree+1)*(1-a)),
                  s*s*hh*common/(1-a),
                  s*hh*jj*common*(degree/(1-a)+a/(1-a)**2)]
        errors = [abs(f-f_n), abs(d-d_n), abs(d2-d2_n)]
        for error, ceiling in zip(errors, bounds):
            assert_bound(error, ceiling)
        tails.append(dict(degree=degree, errors=list(map(float, errors)), bounds=bounds))
    increment = .02*h
    hs, ha = norms(increment)
    a_star, s_star = a+ha, s+hs
    assert a_star < 1
    taylor_error = abs(value(k+increment)-f-first(k, increment))
    taylor_bound = .5*s_star*(2-a_star)*hs**2/(1-a_star)**2
    assert_bound(taylor_error, taylor_bound)
    return dict(label=label, dimension=len(k), schatten_norm=s, operator_norm=a,
                nonnormality=float(np.linalg.norm(k@k.conj().T-k.conj().T@k)),
                commutator=commutator, first_absolute=float(abs(d)), first_bound=bound_d,
                mixed_absolute=float(abs(d2)), mixed_bound=bound_d2,
                symmetry_error=float(symmetry_error), differences=differences, tails=tails,
                taylor_error=float(taylor_error), taylor_bound=float(taylor_bound))


def random_matrices():
    rows = []
    for dim in [3, 7, 12]:
        for target_op in [.15, .5, .8]:
            k = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
            k *= target_op/norms(k)[1]
            h = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
            j = RNG.normal(size=(dim, dim))+1j*RNG.normal(size=(dim, dim))
            h *= .25/norms(h)[1]
            j *= .2/norms(j)[1]
            rows.append(matrix_case(k, h, j, f'complex_nonnormal_{dim}_{target_op}'))
    return rows


def source_cases():
    rows = []
    q, beta, amplitude = 2, 2.3, .25
    lengths = prior.T+np.array([.2, 3.1])**1.5
    refinements = []
    for cutoff in [1, 2, 3]:
        k, v, basis = prior.dirac_matrix(q, cutoff, lengths, beta, amplitude)
        # Color rotations about J_z commute with D0 but rotate the local
        # Hermitian source. In q=2 the retained color basis spans u(2).
        labels = np.array([b[0] for b in basis])
        u = np.diag(np.repeat(np.exp(1j*.7*labels), 2))
        w = np.diag(np.repeat(np.exp(1j*1.9*labels), 2))
        h, j = u@k@u.conj().T/2, w@k@w.conj().T/3
        row = matrix_case(k, h, j, f'fixed_Dirac_q2_cutoff_{cutoff}')
        weighted, radius, integral = prior.integrate_bounds(q, lengths, beta, amplitude, cutoff)
        cq = 4*np.sqrt(2)*np.mean(np.minimum(beta, lengths))/(np.pi*np.sqrt(q))
        # Keep the declared c_q bound; a smaller prior finite-temperature
        # alternative may also be recorded, but is not silently substituted.
        declared_radius = cq*integral
        sk, _ = norms(k)
        assert sk**3 <= declared_radius*(1+5e-10)
        assert declared_radius < 1
        aa = declared_radius**(1/3)
        l3 = integral**(1/3)
        current_ceiling = cq*l3*l3/(1-aa)
        susceptibility_ceiling = cq*l3*(2-aa)/(1-aa)**2
        assert_bound(abs(first(k, h)), current_ceiling*l3/2)
        assert_bound(abs(second(k, j, h)), susceptibility_ceiling*l3*l3/6)
        row.update(lengths=list(map(float, lengths)), beta=beta, amplitude=amplitude,
                   coefficient_cq=float(cq), source_integral_L3_cubed=integral,
                   source_radius=float(declared_radius), previous_smaller_radius=radius,
                   nonlinear_current_each_real_imaginary_norm_bound=float(current_ceiling),
                   nonlinear_current_complex_vector_norm_bound=float(2*current_ceiling),
                   susceptibility_bilinear_norm_bound=float(susceptibility_ceiling))
        if cutoff == 1:
            direct_v = prior.independent_matrix(q, lengths, beta, amplitude, basis)
            error = float(np.linalg.norm(v-direct_v))
            assert error < TOL['rounding']
            # Assemble D0^-1 independently from the basis; compare the
            # local-spacetime quadrature realization with the analytic one.
            inv_d = block_diag(*[-1j*(prior.GAMMA[0]*((2*n+1)*np.pi/beta)+prior.GAMMA[1]*p)
                               /(((2*n+1)*np.pi/beta)**2+p*p) for _, n, p in basis])
            independent_k = inv_d@direct_v
            row['spacetime_multiplier_error'] = error
            row['independent_response_error'] = float(abs(first(k, h)-first(independent_k, h)))
            assert row['independent_response_error'] < TOL['rounding']
        refinements.append(dict(cutoff=cutoff, value_real=float(value(k).real),
                                derivative_real=float(first(k, h).real),
                                hessian_real=float(second(k, j, h).real)))
        rows.append(row)
    return dict(cases=rows, refinements=refinements,
                finite_refinements_prove_continuum_convergence=False,
                continuum_convergence_basis='Ordinary S3 convergence and uniform response tails in note')


def scalar_and_standin():
    counterexamples = []
    for speed in [1, 10, 100, 1000, 10000]:
        k, h = np.array([[.25]]), np.array([[speed/16]])
        derivative = first(k, h).real
        assert abs(derivative-speed/320) < 1e-12
        counterexamples.append(dict(speed=speed, uniform_operator_norm_ceiling=5/16,
                                    derivative=float(derivative)))
    z, v, m = 0., .01, .5
    scalar = m*(1+z)/(1+v*(1+z))
    scalar_prime = m/(1+v*(1+z))**2
    pauli = np.array([[0., 1.], [1., 0.]])
    k, h = scalar*pauli, scalar_prime*pauli
    exact = -2*scalar**3*scalar_prime/(1-scalar**2)
    nonlinear = first(k, h).real
    quadratic = -2*scalar*scalar_prime
    full_relative = -2*scalar*scalar_prime/(1-scalar**2)
    free = 2*v/(1+v)
    full = free+full_relative
    assert abs(nonlinear-exact) < 1e-14
    assert abs(nonlinear+quadratic-full_relative) < 1e-14
    assert abs(free-2/101) < 1e-14
    assert abs(full-(-4798/7701)) < 1e-14
    # Outside the small-op ball: local inverse-gap formula still applies.
    outside_k = np.array([[2.0]])
    outside_h = np.array([[.7]])
    outside_exact = first(outside_k, outside_h).real
    inverse_ceiling = norms(outside_k)[0]**2*np.linalg.norm(np.linalg.inv(np.eye(1)+outside_k), 2)*norms(outside_h)[0]
    assert_bound(abs(outside_exact), inverse_ceiling)
    return dict(value_envelope_counterexample=counterexamples,
                original_standin=dict(scalar=scalar, scalar_prime=scalar_prime,
                                      nonlinear_derivative=float(nonlinear), quadratic_derivative=quadratic,
                                      full_relative_derivative=full_relative, free_derivative=free,
                                      full_derivative=full),
                local_invertible_example=dict(operator_norm=2., derivative=float(outside_exact),
                                              inverse_bound=float(inverse_ceiling)))


def main():
    inputs = ['checks/determinant_response_control.py', 'notes/determinant-response-control.md',
              'checks/nonlinear_current_determinant.py', 'notes/nonlinear-current-determinant.md',
              'checks/source_clock_named_inputs.py', 'notes/source-clock-named-inputs.md',
              'notes/source-clock-unmatched-inputs.md', 'notes/determinant-tail-control.md']
    result = dict(story='S132', date='2026-10-15',
                  status='passed_fixed_domain_nonlinear_current_and_susceptibility',
                  tolerances=TOL, finite_difference_steps=STEPS, series_degrees=DEGREES,
                  seed=20261015, random_nonnormal=random_matrices(),
                  unequal_return_source=source_cases(), scalar_checks=scalar_and_standin(),
                  claims=dict(physical_clock_operator_map_derived=False,
                              interacting_physical_determinant_computed=False,
                              gauge_weighted_current_controlled=False,
                              photon_cutoff_removed=False, lorentzian_continuation=False,
                              s074_gate_closed=False, opening_achieved=False, lean_invoked=False),
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  execution=dict(executable=sys.executable, script=str(Path(__file__).resolve()), cwd=str(Path.cwd())),
                  hash_policy='SHA-256 of file bytes after CRLF-to-LF normalization; no other changes',
                  source_review=dict(url='https://arxiv.org/html/2007.12834v3',
                                     scope='Introduction equation (1.3) and adjacent general complex Schatten scope',
                                     convention='A=-K, k=3', product_anomaly_theorem_imported=False,
                                     novelty_claim=False),
                  input_sha256={name: sha(name) for name in inputs})
    cases = result['random_nonnormal']+result['unequal_return_source']['cases']
    result['summary'] = dict(matrix_cases=len(cases), tail_checks=3*len(DEGREES)*len(cases),
                             max_first_difference_error=max(d['first'] for r in cases for d in r['differences']),
                             max_hessian_difference_error=max(max(d['mixed'], d['path_second']) for r in cases for d in r['differences']),
                             max_hessian_symmetry_error=max(r['symmetry_error'] for r in cases))
    output = ROOT/'evidence/determinant-response-control.json'
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], **result['summary'],
                          execution=result['execution'], hash_policy=result['hash_policy'],
                          script_sha256=result['input_sha256']['checks/determinant_response_control.py'],
                          note_sha256=result['input_sha256']['notes/determinant-response-control.md']), indent=2))


if __name__ == '__main__':
    main()
