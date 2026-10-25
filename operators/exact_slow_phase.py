"""Implicit slow-root boundary primitive; not a global preconditioner.

The root interval encloses an exact, uniquely defined root for each admitted
parameter value. It is never used as a replacement point exponent. The
ordinary domain argument is in notes/current-implicit-phase-graph-domain.md.
Requires python-flint; no finite angular or thermal spectrum is constructed.
"""
from math import comb
from flint import arb


class ExactSlowPhase:
    def __init__(self, L, kappa, K, B0):
        self.L, self.kappa, self.K, self.B0 = map(arb, (L, kappa, K, B0))
        if not (self.L > arb('3e34') and self.kappa > arb('7.8')
                and self.kappa < arb('10.8') and self.K > 0):
            raise ValueError('Outside the proved implicit-root parameter family')
        self.eta = 6/self.L
        low = self.eta*(1-300/self.L)
        high = self.eta*(1+300/self.L)
        # Preserve 6-L*e analytically at the two correlated endpoints.
        rlow = 1800/self.L+self.regular_rest(low)
        rhigh = -1800/self.L+self.regular_rest(high)
        if not (rlow > 0 and rhigh < 0 and -self.L+120 < 0):
            raise ValueError('Implicit-root existence/uniqueness guards failed')
        self.root_bracket = arb((low.lower()+high.upper())/2,
                                (high.upper()-low.lower())/2)
        self.endpoint_signs = (rlow, rhigh)
        self.derivative_interval = -self.L+arb(0, 120)

    def h_regular(self, e):
        e = arb(e)
        if not (e >= 0 and e < arb('1e-30')):
            raise ValueError('Outside proved cotangent remainder interval')
        return ((3+e).digamma()+self.kappa-arb.pi()**2*e/12
                +arb(0, (e**3).abs_upper()))

    def regular_rest(self, e):
        return 5*e+e*e+(6+5*e+e*e)*e*self.h_regular(e)

    def leading_approximation_residual(self):
        """B(4+eta), with eta=6/L. This is NOT set to zero."""
        e = self.eta
        return -self.K*(4+e)*(1+e)*self.regular_rest(e)

    def root_minus_eta_enclosure(self):
        """Mean-value bound for the actual root's offset from eta."""
        r = self.regular_rest(self.eta)
        lo, hi = r/(self.L+120), r/(self.L-120)
        return arb((lo.lower()+hi.upper())/2,(hi.upper()-lo.lower())/2)

    def scaled_pair_jets(self, log_z, order=8):
        """Enclose z^-4 D^j(nu,b), with the EXACT implicit root.

        b=B0(z^4-z^(4+e)), nu=B0(a(4+e)z^(4+e)-a(4)z^4).
        A cutoff equal to one at the end is external. All cutoff derivatives
        and finite-domain returns remain to be included in a global operator.
        """
        x=arb(log_z);e=self.root_bracket
        if not x < 0:
            raise ValueError('The endpoint representation requires 0<z<1')
        p=self.B0*e;exp=(e*x).exp();em1=(e*x).expm1()
        rows=[]
        for j in range(order+1):
            divided=sum((comb(j,k)*arb(4)**(j-k)*e**(k-1)
                         for k in range(1,j+1)),arb(0))
            b=p*(-arb(4)**j*em1/e-divided*exp)
            nu=-arb(7)/15*b-2*p*(4+e)**j*exp/(15*(5+e))
            rows.append((nu,b))
        return rows

    def leading_normal_rows(self):
        """Exact algebraic identities, justified by the certified root.

        The original undifferentiated shape finite-part row is NOT supplied
        by B(e)=0. Returning None prevents treating it as a computed zero.
        """
        return {'trace_pair': arb(0),
                'differentiated_shape_pair': -24*self.K*self.B0,
                'independent_exact_slow_differentiated_shape': arb(0),
                'original_shape_finite_part': None,
                'cutoff_and_global_return': None}

    def scaled_independent_slow_jets(self, log_z, order=8):
        """Enclose the independent exact slow column, before its cutoff.

        This is a separate phase input, not controlled by a fast-input bound.
        Returns z^-4 D^j(nu,b) for (nu,b)=(-a(4+e),1) z^(4+e).
        """
        x=arb(log_z);e=self.root_bracket
        if not x < 0:
            raise ValueError('The endpoint representation requires 0<z<1')
        a=(7+e)/(3*(5+e));exp=(e*x).exp()
        return [(-a*(4+e)**j*exp,(4+e)**j*exp)
                for j in range(order+1)]
