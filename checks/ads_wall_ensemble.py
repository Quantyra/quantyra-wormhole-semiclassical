"""S070: leading fixed-J probe action and boson exchange bookkeeping.

No gravitational determinant, quantum measure or decay rate is computed.
"""
from collections import Counter
from itertools import combinations_with_replacement, permutations
from pathlib import Path
import hashlib
import json
import math
import platform
import mpmath as mp
import sympy as s


def cycle_lengths(perm):
    unseen = set(range(len(perm)))
    lengths = []
    while unseen:
        j = next(iter(unseen))
        length = 0
        while j in unseen:
            unseen.remove(j)
            length += 1
            j = perm[j]
        lengths.append(length)
    return lengths


def symbolic():
    f, J, r, v = s.symbols('f J r v', positive=True)
    p = s.symbols('p', real=True)
    H = s.sqrt(f*f*p*p + f*J*J/(r*r))
    Q = s.sqrt(f + v*v/f)
    ps = s.I*J*v/(r*f*Q)
    # Choose the positive analytic energy branch, continued from v=0.
    Hs = J*f/(r*Q)
    assert s.simplify(Hs**2-H.subs(p, ps)**2) == 0
    assert s.simplify(f*f*ps/Hs-s.I*v) == 0
    assert s.simplify(Hs-s.I*ps*v-J*Q/r) == 0
    ap = s.simplify(s.diff(H, p, 2).subs(p, 0))
    assert s.simplify(ap-r*f**s.Rational(3, 2)/J) == 0
    return {'Euclidean_classical_symbol_action': 'J/r * sqrt(f + r_prime**2/f)',
            'inverse_kinetic_mass': str(ap)}


def trace_checks():
    mp.mp.dps = 60
    rows = []
    for d in (1, 3):
        for beta in (mp.mpf('.3'), mp.mpf('2.1')):
            weights = [mp.exp(-beta*e) for e in (mp.mpf('.4'), mp.mpf('1.3'), mp.mpf('2.7')) for _ in range(d)]
            traces = {c: sum(w**c for w in weights) for c in range(1, 7)}
            for N in range(1, 7):
                fock = sum(mp.fprod(weights[i] for i in inds)
                           for inds in combinations_with_replacement(range(len(weights)), N))
                exchange = sum(mp.fprod(traces[c] for c in cycle_lengths(perm))
                               for perm in permutations(range(N))) / math.factorial(N)
                err = abs(exchange/fock-1)
                assert err < mp.mpf('1e-55')
                rows.append({'N': N, 'degeneracy': d, 'beta': float(beta), 'relative_error': float(err)})
    return rows


def action_checks():
    rows = []
    # Nonlinear fixed-J action, independently differentiated by finite differences.
    for x in (mp.mpf('0'), mp.mpf('.1')):
        R = 2*(1+3*x)
        def f(r):
            return 1+x*r*r-2*(1+2*x)/r+(1+3*x)/(r*r)
        fr = f(R)
        gamma2 = fr*fr/(R*R)-fr*mp.diff(f, R, 2)/2
        mass = 1/(R*fr**mp.mpf('1.5'))
        gamma = mp.sqrt(gamma2)
        for ratio in (mp.mpf('.4'), mp.mpf('1.3')):
            period = 2*mp.pi*ratio/gamma
            for n in (0, 1, 2):
                omega = 2*mp.pi*n/period
                def action(eps):
                    def integrand(t):
                        rr = R+eps*mp.cos(omega*t)
                        vel = -eps*omega*mp.sin(omega*t)
                        return mp.sqrt(f(rr)+vel*vel/f(rr))/rr
                    return mp.quad(integrand, [period*i/16 for i in range(17)])
                norm = period if n == 0 else period/2
                predicted = mass*(omega*omega-gamma2)*norm
                errs = []
                a0 = action(0)
                for h in (mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')):
                    observed = (action(h)+action(-h)-2*a0)/(h*h)
                    errs.append(abs(observed/predicted-1))
                assert errs[-1] < mp.mpf('1e-7')
                assert errs[-1] < errs[0]/9000
                rows.append({'a2_over_L2': float(x), 'period_gamma_over_2pi': float(ratio),
                             'n': n, 'predicted_second_variation': float(predicted),
                             'relative_errors': [float(e) for e in errs]})
    return rows


def mode_checks():
    # Exact rational thresholds; zero eigenvalues must not count as negative.
    checks = 0
    for denominator in (1, 2, 5, 13):
        for numerator in range(1, 26):
            x = s.Rational(numerator, denominator)
            count = sum(1 for n in range(-26, 27) if n*n < x*x)
            assert count == 2*s.ceiling(x)-1
            checks += 1
    sectors = []
    for x in (s.Rational(1, 5), s.Rational(2, 5), s.Integer(1), s.Rational(6, 5)):
        multiplicities = Counter(tuple(sorted(cycle_lengths(p))) for p in permutations(range(4)))
        for lengths, multiplicity in sorted(multiplicities.items()):
            index = sum(int(2*s.ceiling(c*x)-1) for c in lengths)
            zeros = sum(2 for c in lengths if (c*x).is_integer)
            sectors.append({'N': 4, 'beta_gamma_over_2pi': str(x), 'cycle_lengths': lengths,
                            'permutations': multiplicity, 'radial_negative_count': index,
                            'radial_zero_count': zeros})
    return {'exact_threshold_checks': checks, 'exchange_sectors': sectors}


if __name__ == '__main__':
    out = {'date': '2026-10-12', 'status': 'conditional_probe_ensemble_verified',
           'python': platform.python_version(), 'sympy': s.__version__, 'mpmath': mp.__version__,
           'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'symbolic': symbolic(), 'canonical_trace_checks': trace_checks(),
           'nonlinear_action_checks': action_checks(), 'mode_checks': mode_checks(),
           'not_computed': ['exact quantum radial operator or measure', 'full gravitational constrained Hessian',
                            'microscopic BFR radiation state', 'microcanonical mode count',
                            'decay contour or formation rate', 'four-dimensional wormhole stability']}
    target = Path(__file__).resolve().parents[1]/'evidence/ads-wall-ensemble.json'
    target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS: 24 Fock/exchange trace comparisons, 36 nonlinear action differences, 100 exact thresholds')
