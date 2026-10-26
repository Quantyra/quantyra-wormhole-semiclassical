"""Exact mixed row after the original fixed York return, in original units."""
import sympy as sp
from .signed_local_core import SignedLocalCore


class YorkMixedReturn:
    def __init__(self):
        self.core=SignedLocalCore()
        self.z=self.core.z
        self.J,self.br=sp.symbols('J_Y b_radial')
        self.xi=sp.symbols('xi0:6')

    def D(self,e):
        c=self.core;z=self.z;Q,Qp=c.Q[:2]
        out=c.D(e)+sp.diff(e,self.J)*z*(self.br+self.xi[0])/Q
        out+=sp.diff(e,self.br)*((4-z*Qp/Q)*self.br-z*Qp/Q*self.xi[0])
        out+=sum(sp.diff(e,self.xi[j])*self.xi[j+1] for j in range(5))
        return sp.cancel(out)

    def calculate(self):
        c=self.core;z=self.z;Q,Qp=c.Q[:2]
        directions=[[self.xi[0]/2+Q*self.J/z],
                    [-(self.br+3*self.xi[0]+Qp*self.J)/2]]
        for row in directions:
            for _ in range(4):row.append(self.D(row[-1]))
        # 8 pi G A_radius^2 times the PHYSICAL action mixed row.
        mixed=0
        for name in c.densities:
            H=c.build(name)['strong']
            mixed+=sum(-z**4*H[0][b][k]*directions[b][k]
                       for b in range(2) for k in range(5))
        fields=[*self.xi[:5],self.br,self.J]
        generic={str(v):sp.cancel(sp.diff(mixed,v)) for v in fields}
        lam=c.C*c.beta/c.B2
        vacuum={c.lam0:lam+11*lam**2/(360*sp.pi)}
        seed=c.seed_substitution()
        exact={name:sp.cancel(value.subs(vacuum).subs(seed))
               for name,value in generic.items()}
        return {'generic':generic,'seed':exact,'direction_jets':directions}

    def valuation(self,e):
        if e==0:return None
        num,den=sp.fraction(sp.cancel(e))
        def low(p):return min(m[0] for m in sp.Poly(p,self.z).monoms())
        return low(num)-low(den)
