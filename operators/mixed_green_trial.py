"""Explicit mixed-input trial; the original quantum inverse is not implemented here.

Parameters are balls enclosing correlated exact seed quantities. This module
never replaces them by floating point midpoints. The analytic trial is defined
by the formulas below; point evaluations enclose that trial, not an interpolant.
"""
from flint import arb, acb


def cutoff(s, radius):
    """C2 quintic, one through radius and zero beyond twice radius."""
    if s.upper() <= radius.lower():
        return arb(1)
    if s.lower() >= 2 * radius.upper():
        return arb(0)
    x = (s - radius) / radius
    lo = arb(max(arb(0), x.lower()))
    hi = arb(min(arb(1), x.upper()))
    chi = lambda t: 1 - 10*t**3 + 15*t**4 - 6*t**5
    return chi(lo).union(chi(hi))


def legendre_bounds(degree, w, gamma):
    """Full positive-integral bounds for Q and -Q2 d_s Q, degree>=1."""
    assert degree >= 1 and w > 1 and gamma > 0
    eta = w.acosh()
    q = arb.pi() * (-(degree+1)*eta).exp() / (2*(1-(-2*eta).exp()).sqrt())
    flux = arb.pi()*gamma*(degree+1)/4 * (-degree*eta).exp()
    return q, flux


def legendre_q(degree, w):
    """Standard real Q_degree(w), w>1, from its positive integral.

    DLMF 14.12.6 with mu=0 and 14.3.10 gives this normalization.
    The infinite integration tail is added explicitly. The logarithm's
    analytic flag is forwarded to Arb's validated quadrature.
    """
    assert degree >= 1 and w > 1
    upper, _ = legendre_bounds(degree, w, arb(1))
    if upper < arb('1e-80'):
        return arb(0).union(arb(upper.upper())), {'method': 'proved_positive_tail_bound'}
    b = (w*w-1).sqrt()
    end = arb(128)
    def integrand(t, analytic):
        return (-(degree+1)*(acb(w)+b*t.cosh()).log(analytic=analytic)).exp()
    value = acb.integral(integrand, 0, end, rel_tol=arb('1e-18'),
                         abs_tol=arb('1e-60'), eval_limit=30000, depth_limit=60)
    assert value.is_finite() and value.imag.contains(0)
    tail = ((degree+1)*(arb(2).log()-b.log()-end)).exp()/(degree+1)
    out = value.real + arb(0).union(arb(tail.upper()))
    assert out.is_finite()
    return out, {'method': 'validated_positive_integral',
                 'infinite_tail_upper': tail.str(30), 'end': '128'}


class MixedGreenTrial:
    def __init__(self, parameters):
        self.p = parameters
        p = parameters
        self.SD = arb('1e-36')
        self.SP = arb('1e-16')
        self.Ah = p['beta']*p['alpha']*p['L']**2/(2*p['d']*p['gamma'])
        self.At = p['beta']*p['alpha']*p['L']**2/(2*p['gamma']*(1-p['d']*p['v0']))

    def at_bolt(self):
        p = self.p
        psi_n = acb(p['nu']+1).digamma().real
        psi_l = acb(p['ell']+1).digamma().real
        return {
            'h_log_s_coefficient': -self.Ah/2,
            'h_finite_part': self.Ah*((2*p['gamma']/p['a']).log()/2-arb.const_euler()-psi_n),
            't': self.At*(psi_n-psi_l),
        }

    def at(self, s):
        p = self.p
        assert s > 0
        ch, ct = cutoff(s, self.SD), cutoff(s, self.SP)
        if ch.is_zero() and ct.is_zero():
            return {'h': arb(0), 't': arb(0), 'evaluation': {'method': 'exact_cutoff_zero'}}
        w = 1+p['a']*s/p['gamma']
        qn, mn = legendre_q(p['nu'], w)
        ql, ml = legendre_q(p['ell'], w)
        return {'h': ch*self.Ah*qn, 't': ct*self.At*(ql-qn),
                'evaluation': {'heavy': mn, 'slow': ml}}
