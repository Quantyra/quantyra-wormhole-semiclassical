"""Check a first-order energy structure of the published throat equations.

All identities are symbolic. The completed squares agree through order alpha,
not at finite alpha. No spectral stability or physical closure is certified.
Source: arXiv:2610.09847v1, equations 7.3-7.8, 7.26 and Appendix D.
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def matrices(rho, ell, screen, frequency_squared):
    """Return independently transcribed Q1 and Va, and leading N, P1, Vb."""
    f, k = 1 + rho**2, ell - 2
    a, log = sp.atan(rho), sp.log(f)
    n = sp.Matrix([[k, -2, 0], [-2*k, ell+4, ell*screen], [0, 1, screen]])
    p1 = sp.diag(0, -rho/f, 0)
    q1 = sp.Matrix([
        [((8*ell-13)*rho**2+4*ell-5+3*rho*((2*ell-3)*rho**2+4*ell-7)*a
          -(3*ell-6)*log)/3,
         -(20*f-8+3*rho*(5*f+4)*a-6*log)/3, 0],
        [-k*(18*f-5+3*rho*(5*f+4)*a-6*log)/3,
         (6*frequency_squared+f*((14*ell+88)*rho**2+9*ell+68)
          +12*rho*f*((ell+6)*f+ell+4)*a-6*f*(ell+4)*log)/(6*f),
         ell*screen*(14*f-5+12*rho*(f+1)*a-6*log)/6],
        [0, (8*f-4+6*rho*(f+1)*a-3*log)/3,
         screen*(8*f-4+6*rho*(f+1)*a-3*log)/3]])
    va = sp.Matrix([
        [-((8*ell-13)*rho**2+4*ell-5+3*rho*((2*ell-3)*rho**2+4*ell-7)*a
           -(3*ell-6)*log)/3,
         (37*f-13+6*rho*(5*f+4)*a-12*log)/6, 0],
        [k*(39*f-13+6*rho*(5*f+4)*a-12*log)/6,
         -((14*ell+88)*f**2-(5*ell+23)*f+6
           +12*rho*f*((ell+6)*f+ell+4)*a-6*f*(ell+4)*log)/(6*f),
         -ell*screen*(31*f-13+24*rho*(f+1)*a-12*log)/12],
        [0, -(29*f-13+24*rho*(f+1)*a-12*log)/12,
         -screen*(8*f-4+6*rho*(f+1)*a-3*log)/3]])
    vb = sp.diag(0, -1/f, 0)
    return n, p1, q1, va, vb


def verify():
    rho, ell, screen, p2, alpha = sp.symbols('rho L s p2 alpha', real=True)
    f, k = 1+rho**2, ell-2
    n, p1, q1, va, vb = matrices(rho, ell, screen, p2)
    results = {}

    def zero(name, expression):
        value = expression.applyfunc(sp.factor) if isinstance(expression, sp.MatrixBase) else sp.factor(expression)
        assert value == (sp.zeros(*value.shape) if isinstance(value, sp.MatrixBase) else 0), (name, value)
        results[name] = 'exact_zero'
        print(name, 'exact_zero', flush=True)

    g = sp.diag(0, -sp.Rational(1, 2)+1/(2*f), 0)
    q0 = -p2*sp.eye(3)-f*n
    zero('friction_removal', f*g.diff(rho)-p1)
    computed_v1 = -q1 + f*p1.diff(rho)/2 - (g*q0-q0*g)/2
    zero('published_potential_conversion', computed_v1-va-p2*vb)

    h0 = sp.diag(k, 1, ell*screen)
    h1 = sp.diag(-k/sp.Integer(6), 0, -ell*screen/6)
    weighted_v1 = h0*va+h1*f*n
    zero('leading_weighted_symmetry', h0*n-n.T*h0)
    zero('first_order_weighted_symmetry', weighted_v1-weighted_v1.T)
    zero('first_order_kinetic_symmetry', h0*vb-vb.T*h0)
    v = sp.Matrix([-2*screen/k, -screen, 1])
    zero('neutral_leading_mode', n*v)
    norm = (v.T*h0*v)[0]
    den = ell*k+screen*(ell+2)
    published_neutral = -screen/(2*f*den)*(
        8*rho**5*sp.atan(rho)+8*rho**4+16*rho**3*sp.atan(rho)
        +18*rho**2+8*rho*sp.atan(rho)+6-ell*(rho**2-1))
    zero('published_neutral_projection', (v.T*h0*va*v)[0]/norm-published_neutral)
    zero('published_neutral_frequency_coefficient',
         (v.T*h0*vb*v)[0]/norm+screen*k/(f*den))

    # The gravitational entry is R_x/R through first order. The electric
    # entry also matters: a gravitational deformation alone leaves a well.
    s1 = sp.diag((rho+f*sp.atan(rho))/2, rho/(2*f), 0)
    a0 = h0*f*n
    a1 = weighted_v1+h0*f*s1.diff(rho)
    zero('neutral_deformed_first_order', (v.T*a1*v)[0])

    t = 12*rho**3*sp.atan(rho)+17*rho**2+24*rho*sp.atan(rho)-6*sp.log(f)+9
    aa = -t/(12*f)
    bb = -(18*rho**3*sp.atan(rho)+22*rho**2+30*rho*sp.atan(rho)-6*sp.log(f)+17)/(12*f)
    cc = aa+rho**2/(4*f)
    b0 = sp.Matrix([[k, -2, 0], [0, 1, screen]])
    b1 = sp.Matrix([[k*aa, -2*bb, 0], [0, cc, screen*aa]])
    j = sp.diag(1, ell)
    zero('leading_gram_identity', a0-f*b0.T*j*b0)
    zero('first_order_gram_identity', a1-f*(b0.T*j*b1+b1.T*j*b0))
    remainder = ((b0+alpha*b1).T*j*(b0+alpha*b1)
                 -b0.T*j*b0-alpha*(b0.T*j*b1+b1.T*j*b0))
    zero('explicit_second_order_square_remainder', remainder-alpha**2*b1.T*j*b1)
    # Exact finite-alpha positivity does not follow from that identity.
    # For this rational example the truncated deformed matrix is indefinite.
    sample = (a0+alpha*a1).subs({rho: 0, ell: 6, screen: 1})
    determinant = sp.factor(sample.det())
    assert determinant.subs(alpha, sp.Rational(1, 100)) < 0
    assert (a1*v).subs({rho: 0, ell: 6, screen: 1}) != sp.zeros(3, 1)
    return {
        'identities': results,
        'positive_weight_domain': 'L>2, s>0, real rho, 0<=alpha<1',
        'H_alpha_diagonal': ['(L-2)*(1-alpha/6)', '1', 'L*s*(1-alpha/6)'],
        'kinetic_diagonal': ['(L-2)*(1-alpha/6)', '1-alpha/(1+rho^2)', 'L*s*(1-alpha/6)'],
        'sample_truncated_deformed_determinant': str(determinant),
        'sample': 'rho=0, L=6, s=1; determinant negative at alpha=1/100',
        'sample_implication': 'Pointwise indefiniteness at second order; no spectral instability claim',
        'square_identity_order': 'exact coefficients of alpha^0 and alpha^1 only',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = {
        'status': 'passed_first_order_energy_structure_checks',
        'source': 'https://arxiv.org/abs/2610.09847v1',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'checks': verify(),
        'novelty': 'unresolved; established matrix S-deformation method',
        'independent_tensor_derivation_of_throat': False,
        'finite_alpha_positivity_proved': False,
        'global_stability_proved': False,
        'quantum_closure_validity_established': False,
        'lean_certified': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
