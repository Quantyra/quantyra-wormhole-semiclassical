"""Actual finite-seed radial geodesic orbits for the comparison form.

This supplies the global orbit geometry, not an evaluated quantum Hessian
or an approximate inverse. Root intervals enclose exact implicit roots.
The physical shape is even under interchange of the two paired halves.
"""
from math import comb
from flint import arb, acb


def poly(coefficients, z):
    out = z*0
    for c in reversed(coefficients):
        out = out*z+c
    return out


def root_enclosure(f, lo, hi, positive_left=True):
    """Bisection, using a separately proved unique monotone branch."""
    lo, hi = arb(lo), arb(hi)
    direction = 1 if positive_left else -1
    if not (direction*f(lo)>0 and direction*f(hi)<0):
        raise ValueError('Uncertified root bracket')
    for _ in range(400):
        mid=(lo+hi)/2
        v=direction*f(mid)
        if v>0:
            lo=mid
        elif v<0:
            hi=mid
        else:
            break
    return arb((lo.lower()+hi.upper())/2,
               (hi.upper()-lo.lower())/2)


class RadialGeodesicReturn:
    @classmethod
    def from_correlated_seed(cls,beta,u,eps,h):
        """Use the exact root/temperature identities of the retained seed.

        h encloses the original implicit temperature root. This factory does
        not choose a new temperature or assign independent rounded bolt data.
        """
        beta,u,eps,h=map(arb,(beta,u,eps,h))
        v=1+h;zh=1/v
        tau=v*h*h*(beta*(1+2/v+3/v**2)+u/v**2+eps/v**4)
        kval=(v+1)*(u+3*beta*(v*v+1))/(2*v**3)-eps*(v-3)/(2*v**5)
        return cls(beta,u,eps,tau,h*kval,zh)

    def __init__(self, beta, u, eps, tau, gamma, z_h):
        """Inputs must satisfy the exact correlated seed identities.

        Prefer from_correlated_seed; interval consistency alone would not
        establish Q(z_h)=0 or Q'(z_h)=-2 gamma for independent parameters.
        """
        self.beta,self.u,self.eps,self.tau,self.gamma,self.z_h=map(
            arb,(beta,u,eps,tau,gamma,z_h))
        b,u,e,t=self.beta,self.u,self.eps,self.tau
        if not (b>arb('.16') and b<arb('.17') and u>0 and u<arb('1e-10')
                and e<0 and abs(e)<arb('1e-28') and t>0 and t<arb('1e-60')
                and self.gamma>0 and 0<self.z_h and self.z_h<1):
            raise ValueError('Outside the proved finite-seed parameter family')
        if not (6*b+u+e-1).contains(0):
            raise ValueError('Missing exact seed parameter relation')
        # Q(z)=beta+u*z^2-a*z^3+b4*z^4-2*eps*z^5+eps*z^6.
        self.a=4*b+2*u+t
        self.b4=3*b+u+e
        self.coefficients=[b,arb(0),u,-self.a,self.b4,-2*e,e]
        # Q(u*x)-beta=u^3 F(x), evaluated without catastrophic subtraction.
        self.scaled=[arb(0),arb(0),arb(1),-self.a,self.b4*u,-2*e*u*u,e*u**3]
        self.scaled_derivative=[k*self.scaled[k] for k in range(1,7)]
        self.peak_x=root_enclosure(lambda x:poly(self.scaled_derivative,x),'.99','1.01')
        self.separatrix_outer_x=root_enclosure(
            lambda x:poly(self.scaled[2:],x),'1.49','1.51')
        self.peak_z=u*self.peak_x
        self.peak_height_scaled=poly(self.scaled,self.peak_x)

    def Q(self,z):
        return poly(self.coefficients,z)

    def trapped_roots(self,chi):
        """j^2=beta+u^3*chi, for 0<chi<F(peak)."""
        chi=arb(chi)
        if not (chi>0 and chi<self.peak_height_scaled):
            raise ValueError('Trapped orbit requires strict interior chi')
        f=lambda x:poly(self.scaled,x)-chi
        xl=root_enclosure(f,0,self.peak_x,positive_left=False)
        xr=root_enclosure(f,self.peak_x,self.separatrix_outer_x)
        return xl,xr

    def crossing_root(self,fraction):
        """Exact outer root for j^2=fraction*beta; away from separatrix."""
        fraction=arb(fraction)
        if not (fraction>0 and fraction<1):
            raise ValueError('Use the dedicated radial endpoint at j=0')
        # A right bracket short of the bolt avoids cancellation in Q(z_h).
        right=(1+self.separatrix_outer_x*self.u)/2
        while not self.Q(right)-fraction*self.beta<0:
            right=(right+self.z_h)/2
            if not right<self.z_h:
                raise ValueError('Root needs a correlated bolt-end treatment')
        return root_enclosure(lambda z:self.Q(z)-fraction*self.beta,0,right)

    @staticmethod
    def integrate(f,a,b,tol='1e-12'):
        ans=acb.integral(f,a,b,abs_tol=arb(tol),rel_tol=arb(tol),
                         eval_limit=200000,depth_limit=1000)
        if not ans.is_finite() or not ans.imag.contains(0):
            raise ArithmeticError('Orbit integral not certified')
        return ans.real

    def crossing_period(self,fraction,tol='1e-12'):
        zt=self.crossing_root(fraction)
        def integrand(t,analytic):
            z=zt*(1-t*t)
            # Exact Q(zt)=j^2 removes the turning-point zero symbolically.
            g=-sum((self.coefficients[k]*sum(
                (z**i*zt**(k-1-i) for i in range(k)),acb(0))
                for k in range(1,7)),acb(0))
            return 1/g.sqrt(analytic=analytic)
        return 4*zt.sqrt()*self.integrate(integrand,0,1,tol)

    def trapped_period(self,chi,tol='1e-12'):
        xl,xr=self.trapped_roots(chi)
        def integrand(theta,analytic):
            x=(xl+xr)/2+(xr-xl)/2*theta.cos()
            # Both exact root identities remove both endpoint zeros.
            h=-sum((self.scaled[k]*sum(
                (x**i*xl**j*xr**(k-2-i-j)
                 for i in range(k-1) for j in range(k-1-i)),acb(0))
                for k in range(2,7)),acb(0))
            return 1/h.sqrt(analytic=analytic)
        return 2/self.u.sqrt()*self.integrate(integrand,0,arb.pi(),tol)

    def radial_period(self,tol='1e-12'):
        """The j=0 base-profile period; original regular bolt retained."""
        zh,g=self.z_h,self.gamma
        # Use Q(zh)=0 and Q'(zh)=-2 gamma exactly, not interval subtraction.
        c={k:(-1)**k*sum((self.coefficients[n]*comb(n,k)*zh**(n-k)
                         for n in range(k,7)),arb(0)) for k in range(2,7)}
        upper=(zh/(2*g)).sqrt().asinh()
        def integrand(v,analytic):
            sh,ch=v.sinh(),v.cosh()
            normalized=1/(ch*ch)+sum((c[k]*(2*g)**(k-2)*sh**(2*k-2)/(ch*ch)
                                      for k in range(2,7)),acb(0))
            return 1/normalized.sqrt(analytic=analytic)
        return 4*self.integrate(integrand,0,upper,tol)

    def trapped_measure_enclosure(self):
        """Trapped physical-half measure divided by the original time period.

        This is a measure bound, NOT an operator-norm error allowance.
        """
        u,b,e=self.u,self.beta,self.eps
        xcap=arb('1.6')
        delta=xcap*abs(self.b4)*u+2*abs(e)*u*u*xcap*xcap+abs(e)*u**3*xcap**3
        am,ap=self.a-delta,self.a+delta
        if not (am>0 and 1/am<xcap):
            raise ValueError('Cubic squeeze does not cover trapped support')
        fmax=4/(27*am*am)
        scale=16*u**2*u.sqrt()/(15*b.sqrt())
        lower=scale/(ap*ap)*(1-u**3*fmax/(3*b))
        upper=scale/(am*am)
        return arb((lower.lower()+upper.upper())/2,
                   (upper.upper()-lower.lower())/2)

    def frequency(self,period,k):
        """All integer k are retained by the exact spectral representation."""
        return 2*arb.pi()*int(k)/period

    @staticmethod
    def angular_finite_part(t,q):
        """Whole angular tower, in the existing common-restoration scheme."""
        t,q=arb(t),arb(q)
        disc=acb(q*q-4*t).sqrt()
        plus=(q+disc)/2
        minus=2*t/(q+disc)
        val=-2*((1+plus).digamma()+(1+minus).digamma())
        if not (val.is_finite() and val.imag.contains(0)):
            raise ArithmeticError('Angular finite part not enclosed')
        return val.real

    def mass_seam_coefficient(self,j_squared,period,k,mass_coefficient=1):
        """Fourier coefficient of the exact cubic jump subtraction.

        Its mean is zero. Cross terms with the remaining profile are NOT zero.
        It is a decomposition on the flow bundle, not a new metric field.
        """
        j_squared=arb(j_squared)
        if not (j_squared>=0 and j_squared<self.beta):
            raise ValueError('This formula is for crossing trajectories')
        if int(k)==0:
            return arb(0)
        jump=12*arb(mass_coefficient)*(self.beta-j_squared)**arb('1.5')
        return jump/(period*self.frequency(period,k)**4)

    def complete_mass_seam_tail(self,omega,C,q,L,mass_coefficient=1):
        """All crossing trajectories and all |frequency|>omega.

        Returns a comparison quadratic-form upper bound divided by the
        original time period. It is NOT kappa_m or a full column bound.
        No trapped trajectory has this cubic crossing jump. Its remaining
        profile is still part of the full operator.
        """
        om,C,q,L,m=map(arb,(omega,C,q,L,mass_coefficient))
        second_guard=sum((abs(self.coefficients[k])*k*(k-1)
                          for k in range(2,7)),arb(0))
        if not (om>=10 and second_guard<12):
            raise ValueError('Uniform period/frequency guard not established')
        chi_upper=-self.angular_finite_part(C*om*om/4,q)-4*L
        if not chi_upper>2:
            raise ValueError('Monotone logarithmic-tail guard failed')
        a=C/(90*arb.pi()**2)
        bound=(arb(36864)/35*a*self.beta**arb('3.5')
               /(arb.pi()*om**3)*(chi_upper/3+arb(8)/9)*m*m)
        return bound,chi_upper
