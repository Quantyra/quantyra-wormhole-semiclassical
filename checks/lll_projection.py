"""Exact finite-flux LLL projection; the Berezin eigenvalues are established.

Checks the kernel formula independently against normalized monopole-mode
matrix elements. This is a free LLL projection, not full QED matching.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import sympy as sp


def form_factor(q, ell):
    if ell >= q:
        return sp.S.Zero
    return sp.factorial(q) * sp.factorial(q - 1) / (
        sp.factorial(q - ell - 1) * sp.factorial(q + ell))


@lru_cache(None)
def matrix_element(q, index, ell):
    """Expectation of P_ell(u) in one normalized zero mode, index=0..q-1.

    z=(1+u)/2 has Beta(q-index,index+1) density. Compute the integral
    directly from polynomial moments, without using the kernel formula.
    """
    z = sp.Symbol('z')
    poly = sp.Poly(sp.legendre(ell, 2*z-1), z)
    return sum(coef * sp.rf(q-index, power) / sp.rf(q+1, power)
               for (power,), coef in poly.terms())


def verify():
    u, eps = sp.symbols('u eps', real=True)
    kernel_cases = matrix_cases = orthogonality_cases = vertex_cases = 0
    for q in range(1, 13):
        for ell in range(q+2):
            beta = form_factor(q, ell)
            kernel = sp.Poly(sp.expand(((1+u)/2)**(q-1)*sp.legendre(ell, u)), u)
            integral = sum(coef*sp.Rational(2, power+1)
                           for (power,), coef in kernel.terms() if power % 2 == 0)
            assert sp.factor(q*integral/2-beta) == 0
            kernel_cases += 1
            entries = [matrix_element(q, a, ell) for a in range(q)]
            assert sp.factor(sum(x*x for x in entries)-q*beta/(2*ell+1)) == 0
            matrix_cases += 1
            if ell >= q:
                assert all(x == 0 for x in entries)
            # Project the metric vertex independently as a polynomial.
            # m=j-index, v_a=-h_a P_ell'(u) for R=1.
            z = sp.Symbol('z')
            y = sp.legendre(ell, u)
            for a in range(q):
                m = sp.Rational(q-1, 2)-a
                vertex = sp.diff(y, u)*(m-q*u/2)+(1-u*u)*sp.diff(y, u, 2)/4
                poly = sp.Poly(sp.expand(vertex.subs(u, 2*z-1)), z)
                mean = sum(coef*sp.rf(q-a, power)/sp.rf(q+1, power)
                           for (power,), coef in poly.terms())
                assert sp.factor(mean-sp.Rational(ell*(ell+1), 4)*entries[a]) == 0
                vertex_cases += 1
        for ell in range(q):
            for other in range(ell):
                assert sum(matrix_element(q, a, ell)*matrix_element(q, a, other)
                           for a in range(q)) == 0
                orthogonality_cases += 1
    # Fixed-multipole asymptotic, using the product expression directly.
    for ell in range(1, 9):
        beta = sp.prod((1-n*eps)/(1+n*eps) for n in range(1, ell+1))
        assert sp.diff(beta, eps).subs(eps, 0) == -ell*(ell+1)
        assert sp.diff(beta, eps, 2).subs(eps, 0)/2 == sp.Rational(ell**2*(ell+1)**2, 2)
    # The proof of the uniform bound is in the note; these are edge checks.
    for q in range(1, 21):
        for ell in range(q+2):
            beta = form_factor(q, ell)
            assert 0 <= beta <= 1
            assert 1-beta <= sp.Rational(ell*(ell+1), q)
    return {'kernel_integral_cases': kernel_cases,
            'independent_mode_norm_cases': matrix_cases,
            'mode_orthogonality_cases': orthogonality_cases,
            'independent_metric_vertex_integrals': vertex_cases,
            'fixed_multipole_expansions': 8,
            'formula': '(q-1)! q! / ((q-1-l)! (q+l)!) for l<q; zero otherwise',
            'metric_vertex': 'a_eff = a - l(l+1) h/(4 R^2)',
            'screening': 's_l = s * beta_l(q)',
            'scope': 'Free LLL current two-point projection; no higher-Landau-level or full vacuum-stress matching'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = {'status': 'passed_finite_flux_LLL_projection_checks',
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'sympy_version': sp.__version__, 'checks': verify(),
            'novelty': 'Berezin form factor known; specific response/ordering application under review',
            'lean_certified': False, 'full_quantum_closure_proved': False}
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
