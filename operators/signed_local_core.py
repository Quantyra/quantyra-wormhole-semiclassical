"""Both original local strong rows, before adding the canonical determinant.

This is an exact differential-operator generator, not a finite-mode solver.
The remaining canonical determinant has only a shape input/output in this
affine chart. Its common-regulated action must still be supplied. No method
in this module treats that operator, a boundary row or a return as zero.
"""
from math import comb
import sympy as sp


class SignedLocalCore:
    """Use D=z*d/dz and the ORIGINAL dz action measure and W factor."""

    def __init__(self):
        self.z = sp.Symbol('z', positive=True)
        self.C, self.B2 = sp.symbols('C B2', positive=True)
        self.lam0, self.cR, self.qflux, self.gR = sp.symbols(
            'lambda0 cR qflux gR', positive=True)
        self.beta, self.eps, self.tau, self.Fl = sp.symbols('beta eps tau F_lambda')
        self.Q = sp.symbols('Q0:9')
        self.fields = [sp.symbols('nu0:9'), sp.symbols('T0:9')]
        self.zero = {s: 0 for row in self.fields for s in row}
        z, C, B2 = self.z, self.C, self.B2
        nu, T = self.fields
        phi = -sp.log(z)+nu[0]
        dphi = -1+nu[1]
        L = -self.Q[2]-2*self.laplacian(T[0])
        M = self.laplacian(phi)
        N = self.Q[0]*dphi**2/z**2
        curvature = L-6*M-6*N
        self.curvature = curvature
        self.phi, self.L, self.M, self.N = phi, L, M, N
        en2, en4, et2, etm2 = sp.exp(2*nu[0]), sp.exp(4*nu[0]), sp.exp(2*T[0]), sp.exp(-2*T[0])
        self.densities = {
            'Einstein': en2*curvature/(2*z**2)+en2*et2/(C*z**2)
                +3*self.lam0*B2*en4*et2/(C*z**4),
            'Maxwell': -sp.pi*self.qflux**2*et2/(self.gR**2*C*B2),
            'R2': 8*sp.pi*self.cR/(C*B2)*(4*et2+4*C*curvature+C**2*etm2*curvature**2),
            'Weyl': (C**2*etm2*(phi*L**2/60+L*M/30-M**2/10+L*N/30
                      -7*M*N/90-7*N**2/180)-C*phi*L/18
                      +(sp.Rational(1,15)-self.qflux**2/3)*phi*et2-C*N/18)/(2*sp.pi*C*B2),
        }
        # This separately labelled comparison is NOT added to the physical sum.
        self.reference_density = self.Fl*(en2*curvature/(2*z**2)
            +en2*et2/(C*z**2)+3*self.beta*en4*et2/z**4
            -(1/C+3*self.beta+self.eps)*et2
            +4*self.eps*z*sp.exp(-nu[0]+2*T[0])
            -3*self.eps*z**2*sp.exp(-2*nu[0]+2*T[0]))
        self.weyl_endpoint = C*self.Q[0]/(60*sp.pi*B2*z)*(
            phi*self.D(etm2*L)-dphi*etm2*L)
        self._cache = {}

    def D(self, expression):
        """Exact total Euler derivative, including the physical seed jets."""
        result = self.z*sp.diff(expression, self.z)
        for j in range(8):
            result += self.z*self.Q[j+1]*sp.diff(expression, self.Q[j])
            for row in self.fields:
                result += row[j+1]*sp.diff(expression, row[j])
        if any(sp.diff(expression, s) != 0 for s in
               [self.Q[8], self.fields[0][8], self.fields[1][8]]):
            raise ValueError('Requested a jet beyond the declared exact order')
        return sp.expand(result)

    def laplacian(self, expression):
        d = self.D(expression)
        return (self.Q[0]*(self.D(d)-d)+self.z*self.Q[1]*d)/self.z**2

    def adjoint_D(self, expression):
        return -self.D(expression)-expression

    def _adjoint_power(self, expression, n):
        for _ in range(n):
            expression = self.adjoint_D(expression)
        return expression

    def build(self, allocation):
        """Return exact bulk gradient, strong Hessian and boundary jet rows.

        Strong coefficients multiply D^k field, k=0,...,4, in H_c itself.
        Multiplication by z^4 gives the original A=z^4 H_c. Boundary rows
        are retained separately. The fixed physical action factor W is not
        absorbed into a newly named admission norm.
        """
        if allocation in self._cache:
            return self._cache[allocation]
        density = self.reference_density if allocation == 'reference' else self.densities[allocation]
        first = [[sp.expand(sp.diff(density, self.fields[a][i]).subs(self.zero))
                  for i in range(3)] for a in range(2)]
        weak = [[[[sp.expand(sp.diff(density, self.fields[a][i], self.fields[b][j]).subs(self.zero))
                   for j in range(3)] for b in range(2)] for i in range(3)] for a in range(2)]
        source = [sp.expand(sum(self._adjoint_power(first[a][i], i) for i in range(3)))
                  for a in range(2)]
        strong = [[[sp.Integer(0) for _ in range(5)] for _ in range(2)] for _ in range(2)]
        for a in range(2):
            for b in range(2):
                for i in range(3):
                    for j in range(3):
                        for ell in range(i+1):
                            strong[a][b][j+ell] += (-1)**ell*comb(i,ell)*self._adjoint_power(weak[a][i][b][j], i-ell)
                strong[a][b] = [sp.expand(c) for c in strong[a][b]]
        # Integral W_i D^i v dz has lower/upper endpoint difference
        # z[(W1-(D+1)W2)v+W2 Dv]. Never discard this during assembly.
        boundary = [[sp.expand(self.z*(first[a][1]+self.adjoint_D(first[a][2]))),
                     sp.expand(self.z*first[a][2])] for a in range(2)]
        out = {'source':source, 'strong':strong, 'weak':weak, 'boundary_source':boundary}
        self._cache[allocation] = out
        return out

    def apply_symbolic(self, allocation, direction):
        """Apply to arbitrary formal field jets; no metric basis is selected."""
        op = self.build(allocation)['strong']
        return [sp.expand(sum(op[a][b][k]*direction[b][k]
                  for b in range(2) for k in range(5))) for a in range(2)]

    def original_A_coefficients(self, allocation):
        return [[[sp.expand(self.z**4*c) for c in row] for row in block]
                for block in self.build(allocation)['strong']]

    def seed_substitution(self):
        z = self.z
        q = (self.beta+z**2/self.C-(4*self.beta+2/self.C+self.tau)*z**3
             +(3*self.beta+1/self.C+self.eps)*z**4-2*self.eps*z**5+self.eps*z**6)
        return {self.Q[j]:sp.diff(q,z,j) for j in range(9)}

    def independent_trace_variation(self):
        """Linearize the retained Euler TRACE identity, independently of build."""
        z,C,B2 = self.z,self.C,self.B2
        nu,T = self.fields
        phi,L,M,N = self.phi,self.L,self.M,self.N
        en2, en4, et2, etm2 = sp.exp(2*nu[0]), sp.exp(4*nu[0]), sp.exp(2*T[0]), sp.exp(-2*T[0])
        U = 2+C*etm2*self.curvature
        div = lambda expr: (self.D(expr)-expr)/z**2
        FE = en2*self.curvature/z**2+2*en2*et2/(C*z**2)+12*self.lam0*B2*en4*et2/(C*z**4)
        FR = -96*sp.pi*self.cR/B2*(self.laplacian(U)-2*div(self.Q[0]*U*(-1+nu[1])))
        FW = (C**2*(etm2*L**2/60+self.laplacian(etm2*L)/30
            +self.laplacian(etm2*(-M/5-7*N/90))
            -2*div(self.Q[0]*etm2*(L/30-7*M/90-7*N/90)*(-1+nu[1])))
            -C*L/18+(sp.Rational(1,15)-self.qflux**2/3)*et2+C*M/9)/(2*sp.pi*C*B2)
        rows = {'Einstein':FE,'Maxwell':sp.Integer(0),'R2':FR,'Weyl':FW}
        return {name:[[sp.expand(sp.diff(F,row[k]).subs(self.zero)) for k in range(5)]
                      for row in self.fields] for name,F in rows.items()}
