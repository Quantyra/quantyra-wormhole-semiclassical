"""S137.3 exterior response identities; no physical beta closure.

Frozen before execution: exact symbolic 4D constraints/work; 60 digits,
1e-35 quadrature/finite-variation tolerance; fourth-order step 1e-10.
"""
import datetime
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'evidence/exterior-mass-lapse-and-cut-work.json'
LOG = ROOT/'evidence/exterior-mass-lapse-and-cut-work.log'
DPS, TOL, FD_STEP = 60, '1e-35', '1e-10'
INPUTS = ['notes/exterior-mass-lapse-and-cut-work.md',
          'checks/exterior_mass_lapse_and_cut_work.py',
          'notes/exact-spherical-em-reduction.md',
          'notes/sourced-spherical-mass-charge-balance.md',
          'notes/physical-clock-remainder.md',
          'notes/rn-ads-finite-cut-matching.md',
          'notes/exterior-lll-delayed-boundary.md',
          'notes/jt-prescribed-source-constraint.md']


def source_hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')
                             .replace(b'\r', b'\n')).hexdigest() for p in INPUTS}


def record(message):
    with LOG.open('a', encoding='utf-8', newline='\n') as stream:
        stream.write(message+'\n')
    print(message, flush=True)


def zero(expression):
    value = sp.simplify(sp.trigsimp(sp.together(expression)))
    if value != 0:
        value = sp.simplify(sp.expand_trig(value))
    assert value == 0, value


def symbolic_einstein():
    t, r, th, ph = sp.symbols('U R theta phi', real=True)
    G, ell, qm = sp.symbols('G L4 r_m', positive=True)
    alpha = sp.Function('alpha')(t, r)
    F = sp.Function('F')(t, r)
    coords = [t, r, th, ph]
    diagonal = [-sp.exp(2*alpha)*F, 1/F, r*r, r*r*sp.sin(th)**2]
    gamma = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                numerator = ((sp.diff(diagonal[k], coords[i]) if k == j else 0)
                             +(sp.diff(diagonal[k], coords[j]) if k == i else 0)
                             -(sp.diff(diagonal[i], coords[k]) if i == j else 0))
                if numerator != 0:
                    gamma[k, i, j] = sp.simplify(numerator/(2*diagonal[k]))

    def connection(k, i, j):
        return gamma.get((k, i, j), sp.S.Zero)

    def ricci(i, j):
        value = sum(sp.diff(connection(k, i, j), coords[k])
                    -sp.diff(connection(k, i, k), coords[j]) for k in range(4))
        value += sum(connection(k, k, n)*connection(n, i, j)
                     -connection(k, j, n)*connection(n, i, k)
                     for k in range(4) for n in range(4))
        return sp.simplify(value)

    ric = {(i, i): ricci(i, i) for i in range(4)}
    ric[0, 1] = ricci(0, 1)
    zero(ric[3, 3]-sp.sin(th)**2*ric[2, 2])
    scalar = sum(ric[i, i]/diagonal[i] for i in range(4))
    mixed_UU = sp.simplify(ric[0, 0]/diagonal[0]-scalar/2)
    mixed_RR = sp.simplify(ric[1, 1]/diagonal[1]-scalar/2)
    mixed_UR = sp.simplify(ric[0, 1]/diagonal[1])
    closed = (r*sp.diff(F, r)+F-1)/r**2
    zero(mixed_UU-closed)
    zero(mixed_RR-closed-2*F*sp.diff(alpha, r)/r)
    zero(mixed_UR+sp.diff(F, t)/r)

    mass = sp.Function('m')(t, r)
    lapse = 1+r*r/ell**2-2*G*mass/r+qm*qm/r**2
    # Remove Lambda=-3/L4^2 and fixed magnetic stress. Sphere integration
    # converts the remaining 4D Einstein tensor by R^2/(2G).
    rho_equation = r*r/(2*G)*(closed.subs(F, lapse).doit()-3/ell**2+qm*qm/r**4)
    pressure_equation = rho_equation+r*lapse*sp.diff(alpha, r)/G
    flux_equation = -r*sp.diff(lapse, t)/(2*G)
    zero(rho_equation+sp.diff(mass, r))
    zero(flux_equation-sp.diff(mass, t))
    zero(pressure_equation+sp.diff(mass, r)-r*lapse*sp.diff(alpha, r)/G)
    # Exact mass definition, not an additional assumption about its source.
    zero(r/(2*G)*(1+r*r/ell**2+qm*qm/r**2-lapse)-mass)
    constant_mass = sp.symbols('M0', real=True)
    eps = sp.symbols('epsilon', real=True)
    time_profile = sp.Function('A')(t)
    radial_profile = sp.Function('chi')(r)
    zero(rho_equation.subs(mass, constant_mass).doit())
    zero(flux_equation.subs(mass, constant_mass).doit())
    witness_pressure = pressure_equation.subs(mass, constant_mass).doit().subs(
        alpha, eps*time_profile*radial_profile).doit()
    zero(witness_pressure-r*lapse.subs(mass, constant_mass)
         *eps*time_profile*sp.diff(radial_profile, r)/G)
    return dict(full_four_dimensional_ricci_components=5,
                exact_mass_lapse_flux='symbolic zero remainders',
                pressure_only_nonclosure='same mass/flux, distinct cut lapse; no physical matter realization')


def symbolic_work_and_bounds():
    rho, pressure, flux, a, b = sp.symbols('rho p j dotR nR', real=True)
    trace = -rho+pressure
    zero(-rho*a+flux*b-trace*a-(flux*b-pressure*a))
    velocity = sp.symbols('v', real=True)
    gamma2 = 1/(1-velocity**2)
    boosted_p = gamma2*(velocity**2*rho+2*velocity*flux+pressure)
    boosted_j = gamma2*(velocity*rho+(1+velocity**2)*flux+velocity*pressure)
    zero(boosted_j-velocity*boosted_p-(flux+velocity*rho))
    delta = sp.symbols('delta', positive=True)
    P = (delta**-3-1)/3+sp.Rational(3, 2)*(delta**-2-1)+3*(delta**-1-1)+sp.log(1/delta)
    zero(sp.diff(P, delta)+(1+delta)**3/delta**4)
    zero(P.subs(delta, 1))
    # Exact parameter envelope for the classical k=2000 comparison.
    # The retained RN candidate has k^3<=xi<=2k^3; hence the bounds below.
    k, e = Fraction(2000), Fraction(1, 100000)
    mu_max = 4/k**10
    lam_max = mu_max*Fraction(7, 50)
    delta_max = lam_max/e
    assert delta_max == Fraction(546875, 10**34)
    # delta>=1/(4*e*k^10); 4*e*k^10<2^100<exp(100), so log(1/delta)<100.
    assert 4*e*k**10 < 2**100
    upper = (Fraction(8, 3)*e**3+Fraction(3, 2)*mu_max*e**2
             +3*mu_max*lam_max*e+mu_max*lam_max**2*(100+Fraction(49, 4)))
    assert upper < Fraction(2667, 10**18)
    return dict(moving_cut_pressure_work='exact orthonormal and boosted-coordinate identities',
                RN_integral_primitive='exact derivative and endpoint',
                candidate_conditional=dict(k=2000, e=str(e), mu_max=str(mu_max),
                                           delta_max=str(delta_max), alpha_over_A_upper_exact=str(upper),
                                           alpha_over_A_bound='2.667e-15',
                                           physical_chirality_envelope_assigned=False))


def number(value):
    return mp.nstr(value, 55)


def rn_integral_checks():
    rows = []
    for xi, delta in [(mp.mpf(1), mp.mpf(1)/10), (mp.mpf(17), mp.mpf(1)/50)]:
        C = 1+6*xi
        J = lambda s: C+4*xi*s+xi*s*s
        direct = mp.quad(lambda s: (1+s)**3/(s**4*J(s)**2), [delta, 1, mp.inf])
        inverse = mp.quad(lambda v: (v+delta)**3*v**3/
                          (delta**3*(C*v*v+4*xi*delta*v+xi*delta**2)**2), [0, delta, 1])
        error = abs(direct-inverse)/abs(direct)
        assert error < mp.mpf(TOL), error
        P = (delta**-3-1)/3+mp.mpf(3)/2*(delta**-2-1)+3*(delta**-1-1)+mp.log(1/delta)
        bound = P/C**2+1/(4*xi**2)
        assert direct < bound
        travel = mp.quad(lambda s: (1+s)**2/(s*s*J(s)), [delta, 1, mp.inf])
        travel_inverse = mp.quad(lambda v: (v+delta)**2/
                                 (delta*(C*v*v+4*xi*delta*v+xi*delta**2)), [0, delta, 1])
        travel_error = abs(travel-travel_inverse)/travel
        assert travel_error < mp.mpf(TOL), travel_error
        travel_bound = (1/delta+2*mp.log(1/delta)-delta+C/xi)/C
        assert travel <= travel_bound
        rows.append(dict(xi=number(xi), delta=number(delta), I_alpha=number(direct),
                         I_alpha_upper=number(bound), independent_relative_error=number(error),
                         T_ext_over_b=number(travel), T_ext_over_b_upper=number(travel_bound),
                         travel_relative_error=number(travel_error)))
    return rows


def chiral_checks():
    plus = lambda v: mp.exp(-v*v)/1000
    minus = lambda v: -mp.exp(-(v-mp.mpf('.3'))**2/2)/1500
    Pplus = lambda v: mp.sqrt(mp.pi)*mp.erf(v)/2000
    Pminus = lambda v: -mp.sqrt(mp.pi/2)*mp.erf((v-mp.mpf('.3'))/mp.sqrt(2))/1500
    H = lambda U, x: plus(U-x)+minus(U+x)
    current = lambda U, x: -plus(U-x)+minus(U+x)
    U0 = mp.mpf(-2)
    boundary = lambda U: -Pplus(U)+Pplus(U0)+Pminus(U)-Pminus(U0)
    closed_mass = lambda U, x: Pplus(U0)-Pminus(U0)-Pplus(U-x)+Pminus(U+x)

    def quadrature_mass(U, x):
        M = mp.quad(lambda v: current(v, 0), [U0, U], method='gauss-legendre')
        return M-mp.quad(lambda v: H(U, v), [x, 0], method='tanh-sinh')

    step = mp.mpf(FD_STEP)
    derivative = lambda fn, x: (fn(x-2*step)-8*fn(x-step)+8*fn(x+step)-fn(x+2*step))/(12*step)
    rows = []
    hstar = mp.mpf(1)/600
    for U in [mp.mpf('-.2'), mp.mpf('.7'), mp.mpf(2)]:
        for x in [mp.mpf('-.3'), mp.mpf('-.8')]:
            measured = quadrature_mass(U, x)
            errors = [abs(measured-closed_mass(U, x)),
                      abs(derivative(lambda u: quadrature_mass(u, x), U)-current(U, x)),
                      abs(derivative(lambda v: quadrature_mass(U, v), x)-H(U, x))]
            assert max(errors) < mp.mpf(TOL), errors
            assert abs(boundary(U)-measured) <= hstar*abs(x)
            assert abs(boundary(U)) <= hstar*abs(U-U0)
            rows.append(dict(U=number(U), x=number(x), mass=number(measured),
                             quadrature_and_fourth_order_variation_errors=list(map(number, errors))))
    cut = lambda U: -mp.mpf('.5')-mp.mpf('.03')*mp.sin(U)
    cut_derivative = lambda U: -mp.mpf('.03')*mp.cos(U)
    cut_flux = lambda U: current(U, cut(U))+H(U, cut(U))*cut_derivative(U)
    exterior_rate = lambda U: current(U, 0)-cut_flux(U)
    start, stop = mp.mpf('-.2'), mp.mpf('1.4')
    work = mp.quad(exterior_rate, [start, stop], method='gauss-legendre')
    exterior = lambda U: boundary(U)-closed_mass(U, cut(U))
    endpoint_work = exterior(stop)-exterior(start)
    work_error = abs(work-endpoint_work)
    assert work_error < mp.mpf(TOL), work_error
    for U in [start, mp.mpf('.7'), stop]:
        error = abs(derivative(lambda u: closed_mass(u, cut(u)), U)-cut_flux(U))
        assert error < mp.mpf(TOL), error
    return dict(profiles={'J_plus': 'exp(-v^2)/1000',
                          'J_minus': '-exp(-(v-0.3)^2/2)/1500',
                          'scope': 'smooth mathematical verification functions; no quantum-state realization asserted'},
                whole_history_chirality_envelope='1/600', cases=rows,
                integrated_moving_cut_work=number(work), work_error=number(work_error))


def main():
    mp.mp.dps = DPS
    hashes = source_hashes()
    execution = dict(executable=sys.executable, argv=sys.orig_argv, cwd=str(Path.cwd()),
                     command=subprocess.list2cmdline(sys.orig_argv))
    result = dict(milestone='E003/S137.3', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  frozen_protocol=dict(decimal_digits=DPS, absolute_or_relative_tolerance=TOL,
                                       fourth_order_difference_step=FD_STEP), execution=execution,
                  inputs_sha256=hashes,
                  input_hash_policy='SHA-256 after CRLF/CR normalization to LF; no other byte changes',
                  primary_sources=['https://arxiv.org/pdf/gr-qc/9710089',
                                   'https://arxiv.org/html/2308.00871v2'],
                  runtime=dict(python=sys.version, platform=platform.platform(),
                               sympy=sp.__version__, mpmath=mp.__version__))
    record('START '+result['started_utc'])
    record('EXECUTION '+json.dumps(execution))
    record('FROZEN PROTOCOL '+json.dumps(result['frozen_protocol']))
    try:
        result['four_dimensional_constraints'] = symbolic_einstein()
        record('Exact four-dimensional mass/lapse/flux and pressure-only witness: PASS')
        result['moving_work_and_conditional_bound'] = symbolic_work_and_bounds()
        record('Exact moving-cut pressure work and sharp RN primitive/candidate envelope: PASS')
        result['RN_independent_integrals'] = rn_integral_checks()
        record('Independent RN lapse/travel integrals and explicit upper bounds: PASS')
        result['leading_LLL_flux_and_variation'] = chiral_checks()
        record('Leading chiral mass, finite variation and integrated cut flux/work: PASS')
        assert hashes == source_hashes(), 'input changed during execution'
        result['status'] = 'pass_scoped_exterior_response_identities_and_conditional_bounds'
        result['claims'] = dict(exact_spherical_constraints_and_cut_work=True,
                               leading_LLL_response_requires_matched_static_reference=True,
                               corrected_quantum_exterior_or_initial_state_supplied=False,
                               physical_chirality_envelope_established=False,
                               finite_cut_beta_or_O_e_squared_remainder_proved=False,
                               full_preparation_rectangle_certified=False,
                               physical_opening_or_stability_claimed=False, Lean_run=False)
        OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8', newline='\n')
        record('INPUT HASHES '+json.dumps(hashes))
        record('TERMINAL EXIT STATUS 0; scoped component checks PASS; beta remains unclosed')
    except Exception:
        result['status'] = 'failed'
        result['failure_traceback'] = traceback.format_exc()
        OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8', newline='\n')
        record(result['failure_traceback'].rstrip())
        record('FAILURE INPUT HASHES '+json.dumps(hashes))
        record('TERMINAL EXIT STATUS 1; failed attempt retained')
        raise


if __name__ == '__main__':
    main()
