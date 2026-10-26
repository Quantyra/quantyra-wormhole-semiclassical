"""Exact complete angular weights, not an unregulated orbit-trace evaluator.

The common-regulator action identity and original domains are specified in
notes/current-angular-resummed-quantum-operator.md. No thermal/spatial
operator application or strong coupled certificate is supplied here.
"""
from flint import arb, acb


class AngularResponse:
    def __init__(self, C, q, log_mu_B):
        self.C, self.q, self.ell = map(arb, (C, q, log_mu_B))
        if not (self.C > 0 and self.q > 0):
            raise ValueError('Positive original C and q required')

    def _coalescent(self, lam, N=12):
        C, q = self.C, self.q
        r = 1+q/2
        d = q*q-4*C*lam
        t = (abs(d)/(4*r*r)).abs_upper()
        if not t < arb('0.01'):
            raise ValueError('Coalescent series used outside certified disk')
        rs = self.ell-r.digamma()
        bs = arb(0)
        for k in range(1, N+1):
            hz = acb(2*k+1).zeta(acb(r)).real
            rs += (d/4)**k*hz
            bs += k*(d/4)**(k-1)*hz
        re = 4*C*t**(N+1)/(1-t)*(1/r+arb(1)/(2*(N+1)))
        be = 4*C*C/r**2*t**N*((N+1-N*t)/(r*(1-t)**2)+arb(1)/(2*(1-t)))
        return (4*C*rs+arb(0,re.abs_upper()),
                4*C*C*bs+arb(0,be.abs_upper()))

    def resolvent_and_diagonal(self, lam):
        """Enclose FP R_s(lambda) and the entire sum d/(lambda+M^2)^2.

        lambda is a spectral argument; no orbit spectrum is selected here.
        """
        lam = arb(lam)
        if not lam >= 0:
            raise ValueError('Nonnegative orbit spectral argument required')
        C, q = self.C, self.q
        disc = q*q-4*C*lam
        r = 1+q/2
        if (abs(disc)/(4*r*r)).abs_upper() < arb('0.01'):
            return self._coalescent(lam)
        root = acb(disc).sqrt()
        a = 2*C*lam/(q+root)
        b = q-a
        R = 2*C*(2*self.ell-(1+a).digamma()-(1+b).digamma())
        B = 2*C*C*((1+a).polygamma(1)-(1+b).polygamma(1))/root
        # The symmetric functions are real on the entire nonnegative axis.
        if not (R.imag.contains(0) and B.imag.contains(0)):
            raise ArithmeticError('Conjugation symmetry was not enclosed')
        return R.real, B.real

    def bubble(self, lam, mu):
        """Entire angular divided difference, including overlapping inputs."""
        lam, mu = arb(lam), arb(mu)
        if not (lam >= 0 and mu >= 0):
            raise ValueError('Nonnegative orbit spectral arguments required')
        if lam == mu:
            return self.resolvent_and_diagonal(lam)[1]
        if not (mu-lam).contains(0):
            Rl = self.resolvent_and_diagonal(lam)[0]
            Rm = self.resolvent_and_diagonal(mu)[0]
            return (Rl-Rm)/(mu-lam)
        # Monotonicity of the complete positive series provides an enclosure
        # without division by an interval containing zero.
        lo = min(lam.lower(),mu.lower())
        hi = max(lam.upper(),mu.upper())
        lower = self.resolvent_and_diagonal(hi)[1].lower()
        upper = self.resolvent_and_diagonal(lo)[1].upper()
        return arb((lower+upper)/2,(upper-lower)/2)

    def global_bubble_upper(self):
        H2 = arb.pi()**2/6-acb(self.q+1).polygamma(1).real
        return 2*self.C**2*H2/self.q

    def common_heat_regrouping(self):
        """Coefficients of the exact shifted COMMON Mellin identity.

        These are coefficients, not independently renormalized contacts or
        physical operator norms. The three zeta terms AND integral are needed.
        """
        C,q = self.C,self.q
        theta0 = -q-arb(1)/3
        theta1 = (q*q/6-arb(1)/30)/C
        return {'angular_residue':2*C,'theta0':theta0,'theta1':theta1,
                'b0':2*C,'b1':2*C+theta0,'b2':C+theta0+theta1,
                'log_mu_orbit':self.ell-C.log()/2}
