"""Finite resolvent-product field; exact stored fractions, outward queries.

The product representation, not independently rounded partial fractions,
defines the field. T uses the original finite geometry and entire return.
"""
from fractions import Fraction
from math import comb
from flint import arb,acb

def real(q):
    return arb(int(q['numerator']))/arb(int(q['denominator']))

def complex_ball(q):
    return acb(real(q['real']),real(q['imag']))

def polyval(c,x):
    out=x*0
    for a in reversed(c):out=out*x+a
    return out

def polynomials():
    P=[[arb(1)],[arb(0),arb(1)]]
    for j in range(1,4):
        row=[arb(0)]*(j+2)
        for k,v in enumerate(P[j]):row[k+1]+=(2*j+1)*v/(j+1)
        for k,v in enumerate(P[j-1]):row[k]-=j*v/(j+1)
        P.append(row)
    return P

class PairedFilterProfile:
    def __init__(self,parameters,bolt_coefficients,field):
        self.p=parameters;self.bs=bolt_coefficients;self.data=field
        self.a=parameters['a'];self.ga=parameters['gamma'];self.L=parameters['L']
        self.muY=4*self.a/3;self.amplitude=real(field['amplitude'])
        self.base=[acb(real(v)) for v in field['base_poles']]
        terms=[(acb(real(field['constant'])),[]),
               (acb(real(field['large_pole_residue'])),[acb(real(field['large_pole']))])]
        for node in field['upper_half_nodes']:
            pole=complex_ball(node['pole']);weight=complex_ball(node['weight'])
            terms.extend([(weight,[pole]),(weight.conjugate(),[pole.conjugate()])])
        # Aggregate common residues in RAM, always derived from the PRODUCTS.
        common=self.base+[acb(self.muY)]
        weights=[acb(0)]*len(common);individual=[]
        for weight,extra in terms:
            poles=common+extra
            for j,mu in enumerate(poles):
                den=acb(1)
                for k,other in enumerate(poles):
                    if k!=j:den*=other-mu
                assert not den.contains(0)
                coefficient=weight/den
                if j<len(common):weights[j]+=coefficient
                else:individual.append((mu,coefficient))
        self.partials=list(zip(common,weights))+individual
        self.P=polynomials()
        # Original (z Q''-4 Q')' in Legendre polynomial coordinates.
        # Q(z) is reconstructed from its exact bolt expansion.
        Q=[arb(0)]*7
        for j,b in enumerate(bolt_coefficients):
            for k in range(j+1):Q[k]+=b*comb(j,k)*self.L**(j-k)*(-1)**k
        Dz=[(j+2)*(j+1)*(j-3)*Q[j+2] for j in range(5)]
        g=self.ga/self.a;z0=self.L+g;power=[arb(0)]*5
        for j,b in enumerate(Dz):
            for k in range(j+1):power[k]+=b*comb(j,k)*z0**(j-k)*(-g)**k
        leg=[arb(0)]*5
        for j in range(4,-1,-1):
            leg[j]=power[j]/self.P[j][-1]
            for k,b in enumerate(self.P[j]):power[k]-=leg[j]*b
        assert all(v.contains(0) for v in power)
        self.Dleg=leg;self.outer_w=1+self.a*self.L/self.ga
        self.outer_primitive=None

    def kernel(self,mu,w,bolt=False):
        degree=((1+4*mu/self.a).sqrt()-1)/2
        if bolt:
            return -(degree+1).digamma(),acb(0),degree
        def standard_Q(n):
            ww=acb(w)
            pref=arb.pi().sqrt()*((n+1).lgamma()-(n+acb(arb(3)/2)).lgamma()-(n+1)*(2*ww).log()).exp()
            return pref*(1/(ww*ww)).hypgeom_2f1((n+1)/2,(n+2)/2,n+acb(arb(3)/2),abc=True)
        q=standard_Q(degree)
        qplus=standard_Q(degree+1)
        flux=(degree+1)*(qplus-w*q)
        assert q.is_finite() and flux.is_finite()
        return q,flux,degree

    def primitive(self,mu,q,flux,w,bolt=False):
        total=acb(0)
        for j,(dj,Pj) in enumerate(zip(self.Dleg,self.P)):
            den=mu/self.a-j*(j+1);assert not den.contains(0)
            if bolt:numerator=acb(-1)
            else:
                dp=[k*Pj[k] for k in range(1,len(Pj))]
                numerator=polyval(Pj,w)*flux-(w*w-1)*polyval(dp,w)*q
            total+=dj*numerator/den
        return total

    def prepare_return(self):
        if self.outer_primitive is not None:return
        total=acb(0)
        for mu,coefficient in self.partials:
            q,flux,_=self.kernel(mu,self.outer_w)
            total+=coefficient*self.primitive(mu,q,flux,self.outer_w)
        self.outer_primitive=total

    def at_outer(self):
        return {k:arb(0) for k in ['xi','nu','u','J','W']}

    def at(self,s,bolt=False):
        if bolt:s=arb(0)
        else:assert s>0 and s<self.L
        self.prepare_return();w=1+self.a*s/self.ga
        rv=fv=hr=pv=acb(0)
        for mu,coefficient in self.partials:
            q,flux,_=self.kernel(mu,w,bolt=bolt)
            rv+=coefficient*q;fv+=(self.muY-mu)*coefficient*q
            hr+=coefficient*flux
            pv+=coefficient*self.primitive(mu,q,flux,w,bolt=bolt)
        W=self.ga/self.a*(self.outer_primitive-pv)
        vals=[]
        for v in [rv,fv,hr,W]:
            assert v.is_finite() and v.imag.contains(0)
            vals.append(v.real*self.amplitude)
        rv,fv,hr,W=vals;x=s/self.ga;z=self.L-s
        # Geometry is evaluated by the same exact finite-Q jets as the proof.
        from current_source_adapted_fractional_return import J,geometry_jets
        xx=J([x,arb(1)]+[arb(0)]*4)
        zz,q2,R,q=geometry_jets(self.p,self.bs,xx)
        RR=R.c[0];Rx=R.d().c[0];M=z*q.d().d().c[0]+4*self.ga*q.d().c[0]
        nu=z**4/(3*self.L**5)*(-z*RR*fv+(z*RR*self.muY-M)*rv-(7*self.ga*RR-z*Rx)*hr+W)
        u=z**4/self.L**5*(z*RR*fv+(2*M/3-z*RR*self.muY)*rv+(5*self.ga*RR-z*Rx)*hr-2*W/3)
        xi=2*nu+4*self.ga*z**4*RR*hr/(3*self.L**5)
        rx=(self.muY*rv-fv)/2 if bolt else hr/q2.c[0]
        Jvalue=-2*z**5*rx/(3*self.ga*self.L**5)
        return dict(xi=xi,nu=nu,u=u,J=Jvalue,W=W,f=fv,r=rv,hr=hr)
