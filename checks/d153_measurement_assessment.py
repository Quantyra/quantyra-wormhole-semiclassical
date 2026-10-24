"""Independent finite-Fock/operator and sampling checks for the S171 decision.

Development numerics only: this is not a Dirac-field simulation, experiment,
rigorous truncation certificate, or Lean theorem acceptance.
"""
import cmath
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-10
LEVELS = 512
STEPS = 12000


def squeezed(r, phase):
    state = [0j] * LEVELS
    state[0] = 1 / math.sqrt(math.cosh(r))
    for j in range(2, LEVELS, 2):
        state[j] = -cmath.exp(1j*phase)*math.tanh(r)*math.sqrt((j-1)/j)*state[j-2]
    return state


def moments(state):
    n = sum(j*abs(z)**2 for j, z in enumerate(state))
    m = sum(math.sqrt((j+1)*(j+2))*state[j].conjugate()*state[j+2]
            for j in range(len(state)-2))
    return n, m


def quadrature_second(state, theta):
    # Apply X_theta to the Fock vector directly; not the covariance formula.
    z = []
    for j in range(len(state)):
        lower = math.sqrt(j+1)*state[j+1] if j+1 < len(state) else 0j
        upper = math.sqrt(j)*state[j-1] if j else 0j
        z.append((cmath.exp(-1j*theta)*lower+cmath.exp(1j*theta)*upper)/math.sqrt(2))
    return sum(abs(x)**2 for x in z)


def reconstruct(v):
    v0, v45, v90 = v
    n = (v0+v90-1)/2
    return n, complex((v0-v90)/2, v45-n-.5)


def integrate_normal_energy(n, m, tau):
    if tau == 0:
        return n+m.real
    # Sum the number-basis <a†a> and <aa> contractions at each time.
    def integrand(t):
        return (n+(m*cmath.exp(-2j*t)).real)*math.exp(-.5*(t/tau)**2)/(math.sqrt(2*math.pi)*tau)
    left, right = -10*tau, 10*tau
    h = (right-left)/STEPS
    total = integrand(left)+integrand(right)
    total += sum((4 if j % 2 else 2)*integrand(left+j*h) for j in range(1, STEPS))
    return total*h/3


def main():
    errors = []
    cases = []
    norm_error = 0
    for r in [.1, .7, 1.5]:
        for phase in [0, .37, 1.1]:
            state = squeezed(r, phase)
            norm_error = max(norm_error, abs(sum(abs(z)**2 for z in state)-1))
            n, m = moments(state)
            analytic_n = math.sinh(r)**2
            analytic_m = -cmath.exp(1j*phase)*math.sinh(r)*math.cosh(r)
            v = [quadrature_second(state, t) for t in [0, math.pi/4, math.pi/2]]
            n_rec, m_rec = reconstruct(v)
            errors += [abs(n-analytic_n), abs(m-analytic_m), abs(n_rec-n), abs(m_rec-m)]
            for tau in [0, .1, 1, 3]:
                numeric = integrate_normal_energy(n, m, tau)
                analytic = analytic_n+math.exp(-2*tau*tau)*analytic_m.real
                errors.append(abs(numeric-analytic))
                for eta in [.4, .8, 1]:
                    measured = [eta*x+(1-eta)*.5 for x in v]
                    no, mo = reconstruct(measured)
                    observed = no+math.exp(-2*tau*tau)*mo.real
                    errors.append(abs(observed-eta*numeric))
                    assert observed*numeric >= -TOL
                    cases.append(dict(r=r, phase=phase, omega_tau=tau, eta=eta,
                        witness=numeric, observed_witness=observed))
    n, m = moments(squeezed(.7, 0))
    control = dict(v_min=quadrature_second(squeezed(.7, 0), 0),
        local_witness=n+m.real, sampled_witness=integrate_normal_energy(n, m, 1),
        positive_total_excitation=n)
    assert control['v_min'] < .5 and control['local_witness'] < 0 and control['sampled_witness'] > 0
    # A one-photon state: W(0)=-2/pi in the complex-amplitude convention;
    # n=1,m=0 gives positive density for every phase and sampler.
    one = [0j]*LEVELS
    one[1] = 1
    on, om = moments(one)
    assert on == 1 and om == 0
    # Explicit product-POVM illustration: rho_E=diag(.3,.7), M0=|0><0|;
    # sigma_D=|0><0| vs |1><1|. A_D=diag(0,1).
    photon = [.3, .7]
    sigmas = [[1, 0], [0, 1]]
    probabilities = [sum(photon[0]*s for s in sigma) for sigma in sigmas]
    targets = [sigma[1] for sigma in sigmas]
    assert probabilities[0] == probabilities[1] and targets[0] != targets[1]
    maximum = max(errors+[norm_error])
    assert maximum < TOL, maximum
    paths = [Path(__file__), ROOT/'notes/d153-measurement-protocol.md',
        ROOT/'evidence/d153-source-custody.json']
    evidence = dict(scope='finite-Fock and sampler development cross-checks, not physical validation',
        cases=len(cases), levels=LEVELS, simpson_intervals=STEPS, tolerance=TOL,
        max_identity_residual=maximum, max_fock_norm_error=norm_error,
        gaussian_integration_extent_in_tau=10, squeezed_control=control,
        one_photon_control=dict(wigner_at_origin=-2/math.pi, witness=on),
        product_povm_control=dict(probabilities=probabilities, targets=targets),
        cases_detail=cases, input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        experimental_data=False, formal_certificate=False, dirac_source_inferred=False)
    dest = ROOT/'evidence/d153-measurement-assessment.json'
    dest.write_text(json.dumps(evidence, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(f'PASS {len(cases)} scoped cases; max residual {maximum:.3g}')
    print(json.dumps(control, indent=2))


if __name__ == '__main__':
    main()
