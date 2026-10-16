"""S132: finite structural input audit; no clock march or coefficient fit."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes/physical-input-identifiability.md'
SOURCES = [
    'notes/source-clock-running-inputs.md',
    'notes/source-clock-controller.md',
    'checks/source_clock_joint.py',
    'notes/finite-local-electromagnetic-input.md',
    'notes/finite-mouth-determinant.md',
    'checks/finite_local_electromagnetic_input.py',
]
NUMERIC_DPS = 60
FD_STEP = '1e-8'
ABS_TOL = '1e-14'


def zero(value):
    if s.simplify(value) != 0:
        raise AssertionError(str(value))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    # Induction is symbolic in arbitrary positive n; not a finite-order claim.
    mu, mu0, b = s.symbols('mu mu0 b', positive=True)
    C, a, zeta, L = s.symbols('C a zeta_s L', real=True)
    n = s.symbols('n', integer=True, positive=True)
    running = C - b*s.log(mu/mu0)
    nth = (-1)**n*b*s.factorial(n-1)/mu**n
    zero(s.diff(running, mu) - nth.subs(n, 1))
    zero(s.diff(nth, mu) - nth.subs(n, n+1))
    zero(s.diff(nth, C))
    zero(running.subs(C, C+a) - running - a)
    zero(mu*s.diff(running, mu) + b)
    reference_nth = b*(-1)**(n-1)*s.factorial(n-1)/mu0**n
    zero(s.diff(running, mu0) - reference_nth.subs(n, 1))
    zero(s.diff(reference_nth, mu0) - reference_nth.subs(n, n+1))
    zero(s.diff(reference_nth, C))
    # Total reference-scale variation includes dC(mu0)/dmu0=-b/mu0.
    zero(s.diff(running, mu0) - s.diff(running, C)*b/mu0)

    A = s.symbols('A', positive=True)
    aggregate = A*(zeta+2*C*L)
    row = s.Matrix([[s.diff(aggregate, zeta), s.diff(aggregate, C)]])
    if row.rank() != 1:
        raise AssertionError('one aggregate must have rank one')
    zero((row*s.Matrix([-2*L, 1]))[0])
    shifted = aggregate.subs({zeta: zeta-2*L*a, C: C+a}, simultaneous=True)
    zero(shifted-aggregate)

    A1, A2 = s.symbols('A1 A2', nonzero=True, real=True)
    L1, L2, q1, q2 = s.symbols('L1 L2 q1 q2', real=True)
    matrix = s.Matrix([[A1, 2*A1*L1], [A2, 2*A2*L2]])
    zero(matrix.det() - 2*A1*A2*(L2-L1))
    inverse_C = (q2-q1)/(2*(L2-L1))
    inverse_zeta = (L2*q1-L1*q2)/(L2-L1)
    zero(inverse_zeta+2*L1*inverse_C-q1)
    zero(inverse_zeta+2*L2*inverse_C-q2)
    if matrix.subs(L2, L1).rank() != 1:
        raise AssertionError('equal weights must remain rank one')

    cJ, B0, kg, chi_data = s.symbols('cJ B0 kappa_grav chi_data', real=True)
    chi = 1+2*cJ*kg**4*B0**2
    zero((chi-1)/(2*kg**4*B0**2)-cJ)
    zero(chi.subs(B0, 0)-1)
    zero(s.diff(chi.subs(cJ, 2), cJ))
    zero(s.diff(chi, B0, 2)-4*cJ*kg**4)
    J = s.symbols('J', real=True)
    local = cJ*kg**4*J**2/16
    zero(local.subs(J, 0))
    zero(s.diff(local, J).subs(J, 0))

    r0, kc, z, g, Dj = s.symbols('r0 kappa_clock z g D_j', real=True)
    acceleration = g/2-(1+aggregate*r0**2)*z-4*aggregate*kc*r0*z**2-3*aggregate*kc**2*z**3
    zero(s.diff(acceleration, Dj))
    zero(acceleration.subs({zeta:zeta-2*L*a, C:C+a}, simultaneous=True)-acceleration)
    eta, k, T = s.symbols('eta k T', positive=True)
    mouth = eta*s.coth(k*(T+Dj)/2)/(2*k)
    zero(s.diff(mouth, Dj)+eta/(4*s.sinh(k*(T+Dj)/2)**2))

    # Independent exact rational forward and inverse calculations.
    # Both couplings differ and all toy values are nonzero.
    pair1, pair2 = (F(13), F(2)), (F(7), F(3))
    observed1 = F(2)*(pair1[0]+2*pair1[1]*F(3))
    observed2 = F(2)*(pair2[0]+2*pair2[1]*F(3))
    if observed1 != observed2 or observed1 != 50:
        raise AssertionError('rational ambiguity witness failed')
    # Known synthetic geometries: A=(2,3), L=(3,5), common zeta=13,C=2.
    lambda1, lambda2 = F(50), F(99)
    datum1, datum2 = lambda1/F(2), lambda2/F(3)
    recovered_C = (datum2-datum1)/(2*(F(5)-F(3)))
    recovered_zeta = (F(5)*datum1-F(3)*datum2)/(F(5)-F(3))
    if (recovered_zeta, recovered_C) != pair1:
        raise AssertionError('rational two-geometry inversion failed')
    toy_cJ, toy_B, toy_kg = F(7,5), F(3,2), F(2,3)
    toy_chi = 1+2*toy_cJ*toy_kg**4*toy_B**2
    recovered_cJ = (toy_chi-1)/(2*toy_kg**4*toy_B**2)
    if recovered_cJ != toy_cJ or toy_chi != F(101,45):
        raise AssertionError('rational susceptibility inversion failed')
    # D_j values affect a full mouth response, but not the retained clock RHS.
    def rational_clock(dj):
        return F(1,2)-(1+observed1*F(1,100)**2)*F(1,4)-4*observed1*F(1,50)*F(1,100)*F(1,4)**2-3*observed1*F(1,50)**2*F(1,4)**3
    if rational_clock(F(1,3)) != rational_clock(F(4,3)):
        raise AssertionError('retained clock length invariance failed')

    mp.mp.dps = NUMERIC_DPS
    h, tolerance = mp.mpf(FD_STEP), mp.mpf(ABS_TOL)
    beta = 1/(20*mp.pi**2)
    x, x0 = mp.mpf('1.5'), mp.mpf('0.5')
    def run(constant, scale):
        return constant-beta*mp.log(scale/x0)
    derivatives = [(run(constant,x+h)-run(constant,x-h))/(2*h) for constant in (mp.mpf(2),mp.mpf(3))]
    run_errors = [abs(value+beta/x) for value in derivatives]
    if max(run_errors) >= tolerance or abs(run(3,x)-run(2,x)-1) >= tolerance:
        raise AssertionError('independent scale check failed')
    def run_reference(constant, reference):
        return constant-beta*mp.log(x/reference)
    reference_derivatives = [(run_reference(constant,x0+h)-run_reference(constant,x0-h))/(2*h) for constant in (mp.mpf(2),mp.mpf(3))]
    reference_errors = [abs(value-beta/x0) for value in reference_derivatives]
    if max(reference_errors) >= tolerance:
        raise AssertionError('independent fixed-C reference-scale partial failed')
    # Vary background, retaining the unknown coefficient in the forward law.
    def susceptibility(background):
        return 1+2*(mp.mpf(7)/5)*(mp.mpf(2)/3)**4*background**2
    center = mp.mpf(3)/2
    fd_chi = (susceptibility(center+h)-susceptibility(center-h))/(2*h)
    susceptibility_error = abs(fd_chi-4*(mp.mpf(7)/5)*(mp.mpf(2)/3)**4*center)
    if susceptibility_error >= tolerance:
        raise AssertionError('independent susceptibility variation failed')
    mouth_values = [mp.coth((1+dj)/2)/2 for dj in (mp.mpf(1)/3,mp.mpf(4)/3)]
    if abs(mouth_values[0]-mouth_values[1]) <= tolerance:
        raise AssertionError('full mouth model must retain length dependence')

    return {
        'verdict':'structural ambiguity established; derivative-only continuation closed',
        'all_order_scale_derivative':str(nth),
        'all_order_fixed_C_mu0_partial':str(reference_nth),
        'all_order_proof':'exact symbolic base and arbitrary-positive-integer induction step',
        'reference_scale_scope':'mu0 partial holds C independent; physical total derivative adds dC/dmu0=-b/mu0 and vanishes',
        'aggregate_rank':1,
        'aggregate_kernel':['-2*L','1'],
        'two_geometry_determinant':str(s.factor(matrix.det())),
        'two_geometry_assumptions':'same scale, scheme and common couplings; independently known nonzero A_i and distinct L_i; independently matched Lambda_i',
        'rational_same_data_witness':{'zeta_C_1':list(map(str,pair1)),'zeta_C_2':list(map(str,pair2)),'A':'2','L':'3','Lambda':str(observed1)},
        'rational_two_geometry_inverse':{'A':['2','3'],'L':['3','5'],'Lambda':[str(lambda1),str(lambda2)],'recovered_zeta':str(recovered_zeta),'recovered_C':str(recovered_C)},
        'rational_susceptibility_inverse':{'cJ':str(toy_cJ),'B0':str(toy_B),'kappa_grav':str(toy_kg),'chi':str(toy_chi),'recovered_cJ':str(recovered_cJ)},
        'numeric_checks':{'dps':NUMERIC_DPS,'step':FD_STEP,'absolute_tolerance':ABS_TOL,'scale_derivative_absolute_errors':list(map(str,run_errors)),'fixed_C_mu0_partial_absolute_errors':list(map(str,reference_errors)),'susceptibility_absolute_error':str(susceptibility_error),'mouth_values':list(map(str,mouth_values))},
        'physical_measurements_supplied':False,
        'synthetic_examples_only':True,
        'S074_complete':False,
        'new_Lambda_interval':False,
        'area_budget_evaluated':False,
        'opening_achieved':False,
        'stability_achieved':False,
        'transport_achieved':False,
        'preparation_achieved':False,
        'creation_achieved':False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'evidence/physical-input-identifiability.json')
    args = parser.parse_args()
    result = verify()
    inputs = SOURCES + ['notes/physical-input-identifiability.md','checks/physical_input_identifiability.py']
    result['inputs_sha256'] = {name:sha(ROOT/name) for name in inputs}
    result['runtime'] = {'python':sys.version,'platform':platform.platform(),'sympy':s.__version__,'mpmath':mp.__version__}
    result['status'] = 'pass'
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('S132 physical input identifiability: PASS')
    print('All-order constant ambiguity, clock rank/kernel, conditional inverse, susceptibility and length checks passed.')
    print('Examples are synthetic; no physical input value or full S074 closure follows.')
    print('Receipt: '+str(args.output))


if __name__ == '__main__':
    main()
