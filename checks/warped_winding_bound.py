"""Cross-check the ordinary warped-circle fermion winding bound.

The proof is in notes/warped-winding-bound.md. Numerical checks support
normalization, local metric variation and transfer evaluation; they do not
certify all floating-point errors or a real MMP exterior.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm, expm_frechet
from scipy.special import k1
import sympy as sp


SIGMA_A = np.diag([-1.0, 1.0])
SIGMA_B = np.array([[0.0, 1.0], [1.0, 0.0]])


def transfer_data(lam, frequency, segments, derivative=None):
    """Segments are (coordinate width, lapse N, radial B, radius R)."""
    matrices = []
    derivatives = []
    log_scale = 0.0
    for index, (width, lapse, radial, radius) in enumerate(segments):
        a, b = lam*radial/radius, frequency*radial/lapse
        generator = a*SIGMA_A+b*SIGMA_B
        k = np.hypot(a, b)
        scale = k*width
        # Remove a scalar exponential for stability. It cancels in derivative ratios.
        shifted = generator*width-scale*np.eye(2)
        if derivative is not None and index == derivative[0]:
            which = derivative[1]
            direction = {'N': -b/lapse*SIGMA_B,
                         'B': a/radial*SIGMA_A+b/radial*SIGMA_B,
                         'R': -a/radius*SIGMA_A}[which]*width
            mat, dmat = expm_frechet(shifted, direction)
        else:
            mat, dmat = expm(shifted), np.zeros((2, 2))
        matrices.append(mat)
        derivatives.append(dmat)
        log_scale += scale
    product, dproduct = np.eye(2), np.zeros((2, 2))
    for mat, dmat in zip(matrices, derivatives):
        dproduct = dmat@product+mat@dproduct
        product = mat@product
    trace = np.trace(product)
    log_trace = np.log(trace)+log_scale
    root = np.sqrt(max(0.0, 1-4*np.exp(-2*log_trace)))
    zeta = log_trace+np.log((1+root)/2)
    dzeta = np.trace(dproduct)/trace/root if derivative is not None else None
    return zeta, dzeta


def energy(lam, segments):
    length = sum(w*b/n for w, n, b, _ in segments)
    value, error = quad(lambda z: -2/np.pi*np.log1p(np.exp(-transfer_data(lam, z, segments)[0])),
                        0, 40/length, epsabs=2e-12, epsrel=2e-11)
    return value, error+2*np.exp(-40)/(np.pi*length)


def metric_derivative(lam, segments, index, which):
    length = sum(w*b/n for w, n, b, _ in segments)
    def integrand(frequency):
        zeta, derivative = transfer_data(lam, frequency, segments, (index, which))
        return 2/np.pi*derivative*np.exp(-zeta)/(1+np.exp(-zeta))
    return quad(integrand, 0, 40/length, epsabs=2e-12, epsrel=2e-10)[0]


def local_kernels(lam, frequency, segments, index):
    """Insert each generator at the midpoint of a segment, using cyclicity."""
    width, lapse, radial, radius = segments[index]
    half = (width/2, lapse, radial, radius)
    cyclic = [half]+segments[index+1:]+segments[:index]+[half]
    product, log_scale = np.eye(2), 0.0
    for width, lapse, radial, radius in cyclic:
        a, b = lam*radial/radius, frequency*radial/lapse
        scale = np.hypot(a, b)*width
        product = expm(width*(a*SIGMA_A+b*SIGMA_B)-scale*np.eye(2))@product
        log_scale += scale
    trace = np.trace(product)
    log_trace = np.log(trace)+log_scale
    root = np.sqrt(max(0.0, 1-4*np.exp(-2*log_trace)))
    zeta = log_trace+np.log((1+root)/2)
    return zeta, np.trace(SIGMA_A@product)/trace/root, np.trace(SIGMA_B@product)/trace/root


def relative_bounds(q, barrier, local_kappa):
    first = np.sqrt(q+1)
    I = barrier*first
    assert I >= 3
    energy_poly = ((q+2)*(1+I)
                   + 2*(first**2+3*first/barrier+3/barrier**2))
    gamma = local_kappa/barrier
    p = (gamma+1)*I**2+(gamma+2)*I+2
    primitive = ((gamma+1)*(I**3+3*I**2+6*I+6)
                 +(gamma+2)*(I**2+2*I+2)+2*(I+1))
    null_poly = (q+2)*p+2*primitive/barrier**2
    return {'energy_ratio_upper_bound': float(24*np.exp(-I)*energy_poly/(np.pi**2*q)),
            'null_stress_ratio_upper_bound': float(12*np.exp(-I)*null_poly/(np.pi**2*q*np.tanh(I)))}


def spectral_check():
    """Independent heat-regulated Hamiltonian spectrum for m(x)=0.7+0.3 cos(x)."""
    length = 2*np.pi
    def log_winding(frequency):
        def ode(x, flat):
            return (((0.7+0.3*np.cos(x))*SIGMA_A+frequency*SIGMA_B)
                    @flat.reshape(2, 2)).ravel()
        result = solve_ivp(ode, [0, length], np.eye(2).ravel(),
                           method='DOP853', rtol=2e-11, atol=2e-12)
        assert result.success
        zeta = np.arccosh(np.trace(result.y[:, -1].reshape(2, 2))/2)
        return -2/np.pi*np.log1p(np.exp(-zeta))
    transfer_energy, quad_error = quad(log_winding, 0, 32/length, epsabs=2e-12, epsrel=2e-10)
    # A covariant energy regulator is essential; a bare Fourier-index cutoff
    # alone would not respect large-gauge equivalence of the phase sectors.
    cutoff = 100
    modes = np.arange(-cutoff, cutoff+1)
    count = len(modes)
    mass = 0.7*np.eye(count)+0.15*(np.eye(count, k=1)+np.eye(count, k=-1))
    alpha = np.array([[0, -1j], [1j, 0]])
    temperatures = np.array([0.04, 0.02, 0.01, 0.005])
    def regulated(phase):
        hamiltonian = np.kron(alpha, np.diag(modes+phase))+np.kron(SIGMA_B, mass)
        levels = np.linalg.eigvalsh(hamiltonian)
        return -0.5*np.sum(np.abs(levels)[:, None]*np.exp(-levels[:, None]**2*temperatures), axis=0)
    phase_points = 24
    average = sum(regulated(-0.5+(j+0.5)/phase_points) for j in range(phase_points))/phase_points
    differences = regulated(0.5)-average
    extrapolations = [(differences[j]-6*differences[j+1]+8*differences[j+2])/3 for j in [0, 1]]
    residuals = [abs(x-transfer_energy) for x in extrapolations]
    assert residuals[1] < 2e-7 and residuals[1] < residuals[0]/4
    return {'profile': 'm(x)=0.7+0.3 cos(x), optical circumference 2 pi',
            'Fourier_cutoff': cutoff, 'phase_quadrature_points': phase_points,
            'heat_regulator_parameters': temperatures.tolist(),
            'regulated_winding_energies': differences.tolist(),
            'quadratic_extrapolations': [float(x) for x in extrapolations],
            'transfer_energy': float(transfer_energy),
            'transfer_quadrature_estimate_plus_tail': float(quad_error+2*np.exp(-32)/(np.pi*length)),
            'spectral_vs_transfer_residuals': [float(x) for x in residuals],
            'spectral_truncation_and_extrapolation_not_certified': True}


def verify():
    symbolic = {}
    z, gamma, T, lam = sp.symbols('z gamma T lambda', positive=True)
    p = (gamma+1)*z**2+(gamma+2)*z+2
    primitive = ((gamma+1)*(z**3+3*z**2+6*z+6)
                 +(gamma+2)*(z**2+2*z+2)+2*(z+1))
    expressions = {
        'stress_tail_primitive': sp.diff(-sp.exp(-z)*primitive, z)-z*p*sp.exp(-z),
        'energy_tail_primitive': sp.diff(-2*sp.exp(-T*lam)*(lam**2+3*lam/T+3/T**2), lam)
                                 -2*lam*(1+T*lam)*sp.exp(-T*lam),
        'zeroth_frequency_moment_derivative': sp.diff((z+1)*sp.exp(-z), z)+z*sp.exp(-z),
        'first_frequency_moment_derivative': sp.diff((z*z/2+z+1)*sp.exp(-z), z)+z*z*sp.exp(-z)/2
    }
    for name, expression in expressions.items():
        assert sp.simplify(expression) == 0, name
        symbolic[name] = 'exact_zero'

    # Constant mass reproduces the independent Bessel-series cylinder result.
    constant = []
    for mass, length in [(0.0, 3.0), (0.3, 2.0), (1.0, 3.0), (2.0, 4.0)]:
        value, error = energy(mass, [(length, 1.0, 1.0, 1.0)])
        expected = (-np.pi/(6*length) if mass == 0 else
                    2*mass/np.pi*sum((-1.0)**k*k1(k*mass*length)/k for k in range(1, 241)))
        assert abs(value-expected) < 3e-11
        constant.append({'mass': mass, 'length': length, 'energy': value,
                         'Bessel_or_massless_residual': abs(value-expected), 'estimated_error_plus_tail': error})

    profiles = [
        [(1.0, 1.0, 1.0, 1.0), (8.0, 1.0, 1.0, 100.0)],
        [(0.7, 0.6, 1.4, 0.9), (1.3, 1.5, 0.8, 3.0), (0.4, 0.9, 1.1, 1.2)],
        [(1.0, 2.0, 0.7, 2.0), (2.0, 0.7, 1.1, 0.8)]
    ]
    trace_cases, kernel_cases = 0, 0
    for profile in profiles:
        length = sum(w*b/n for w, n, b, _ in profile)
        barrier = sum(w*b/r for w, _, b, r in profile)
        for mass in [0.3, 1.0, 2.0, 4.0]:
            for scaled_frequency in [0.0, 0.1, 1.0, 4.0, 12.0]:
                zeta, _ = transfer_data(mass, scaled_frequency/length, profile)
                assert zeta+2e-12 >= max(mass*barrier, scaled_frequency)
                trace_cases += 1
                for index in range(len(profile)):
                    zeta_local, ka, kb = local_kernels(mass, scaled_frequency/length, profile, index)
                    assert abs(zeta_local-zeta) < 3e-12
                    assert min(ka, kb) >= -2e-12
                    assert max(ka, kb) <= 1/np.tanh(zeta)+2e-12
                    kernel_cases += 1
        value, _ = energy(2.0, profile)
        I = 2*barrier
        assert 0 < -value <= 2*(I+1)*np.exp(-I)/(np.pi*length)

    # Differentiate the determinant, then compare to independently varied metrics.
    variations = []
    for profile in profiles[:2]:
        for index in range(len(profile)):
            for which, position in [('N', 1), ('B', 2), ('R', 3)]:
                analytic = metric_derivative(2.0, profile, index, which)
                step = 2e-4*profile[index][position]
                plus, minus = [list(row) for row in profile], [list(row) for row in profile]
                plus[index][position] += step
                minus[index][position] -= step
                finite = (energy(2.0, plus)[0]-energy(2.0, minus)[0])/(2*step)
                residual = abs(analytic-finite)
                assert residual < 2e-8
                expected_sign = 1 if which == 'B' else -1
                assert expected_sign*analytic > 0
                variations.append({'segment': index, 'parameter': which, 'derivative': analytic,
                                   'finite_difference_residual': residual})

    # Smooth AdS throat: independent ODE compared with a refined product.
    rho_cutoff, mass, frequency = 2.0, 2.0, 0.7
    xc = np.arctan(rho_cutoff)
    def ode(x, flat):
        return ((mass/np.cos(x)*SIGMA_A+frequency*SIGMA_B)@flat.reshape(2, 2)).ravel()
    ode_result = solve_ivp(ode, [-xc, xc], np.eye(2).ravel(), rtol=2e-11, atol=2e-12)
    assert ode_result.success
    zeta_ode = np.arccosh(np.trace(ode_result.y[:, -1].reshape(2, 2))/2)
    refinements = []
    for count in [128, 256, 512]:
        width = 2*xc/count
        smooth = [(width, 1/np.cos(-xc+(j+0.5)*width),
                   1/np.cos(-xc+(j+0.5)*width), 1.0) for j in range(count)]
        measured, _ = transfer_data(mass, frequency, smooth)
        refinements.append(abs(measured-zeta_ode))
    assert refinements[2] < refinements[1]/3 < refinements[0]/9
    assert zeta_ode >= max(2*mass*np.arcsinh(rho_cutoff), 2*xc*frequency)

    # The tail bounds come from the ordinary monotonicity/integral proof, not a truncation.
    samples = []
    for q in [3, 8, 16, 64]:
        for rho in [2.0, 8.0]:
            barrier = 2*np.arcsinh(rho)
            bounds = relative_bounds(q, barrier, np.pi)
            I = barrier*np.sqrt(q+1)
            gamma_num = np.pi/barrier
            summed_e = summed_s = 0.0
            for n in range(1, 201):
                v = barrier*np.sqrt(n*(n+q))
                summed_e += (q+2*n)*(1+v)*np.exp(-v)
                summed_s += (q+2*n)*((gamma_num+1)*v*v+(gamma_num+2)*v+2)*np.exp(-v)/np.tanh(v)
            assert 24*summed_e/(np.pi**2*q) <= bounds['energy_ratio_upper_bound']
            assert 12*summed_s/(np.pi**2*q) <= bounds['null_stress_ratio_upper_bound']
            samples.append({'q': q, 'rho_cutoff': rho, 'integrated_barrier_lower_bound': barrier,
                            'local_NL_over_R': float(np.pi), **bounds})
    return {'symbolic': symbolic, 'constant_profile_checks': constant,
            'transfer_trace_bound_cases': trace_cases, 'local_kernel_bound_cases': kernel_cases,
            'local_metric_variations': variations,
            'smooth_AdS_transfer_refinement_errors': refinements, 'all_mode_bound_samples': samples,
            'independent_spectral_winding_check': spectral_check()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'status': 'passed_warped_winding_cross_checks',
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'checks': verify(),
              'ordinary_energy_and_local_null_stress_bound': True,
              'reference': 'Antiperiodic ground-state energy/stress minus the average over spatial holonomy on the identical prescribed warped circle',
              'angular_modes_assumed_unmixed': True, 'gauge_field_prescribed': True,
              'full_numerical_error_certified': False,
              'actual_MMP_exterior_matched': False, 'full_support_stress_controlled': False,
              'physical_stability_proved': False, 'novelty_established': False,
              'lean_certified': False, 'lean_invocations': 0}
    args.output.write_bytes((json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps({'status': result['status'], 'script_sha256': result['script_sha256'],
                      'transfer_cases': result['checks']['transfer_trace_bound_cases'],
                      'metric_variations': len(result['checks']['local_metric_variations']),
                      'bound_samples': result['checks']['all_mode_bound_samples']}))


if __name__ == '__main__':
    main()
