"""One-Dirac-loop trace matching for magnetic AdS2 x S2 throats.

Independent metric curvature and Clifford/heat-kernel contractions. This
tests a necessary background equation, not a complete stress or solution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def curvature(metric,coords):
    inverse=metric.inv()
    Gamma={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                Gamma[a,b,c]=sp.simplify(sum(inverse[a,d]*(sp.diff(metric[d,c],coords[b])
                    +sp.diff(metric[d,b],coords[c])-sp.diff(metric[b,c],coords[d]))/2 for d in range(4)))
    riemann={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    riemann[a,b,c,d]=sp.simplify(sp.diff(Gamma[a,d,b],coords[c])-sp.diff(Gamma[a,c,b],coords[d])
                        +sum(Gamma[a,c,e]*Gamma[e,d,b]-Gamma[a,d,e]*Gamma[e,c,b] for e in range(4)))
    Ric=sp.Matrix(4,4,lambda b,d:sp.simplify(sum(riemann[a,b,a,d] for a in range(4))))
    scalar=sp.simplify(sp.trace(inverse*Ric))
    ric2=sp.simplify(sp.trace(inverse*Ric*inverse*Ric))
    assert metric.is_diagonal()
    riem2=sp.simplify(sum(metric[a,a]*inverse[b,b]*inverse[c,c]*inverse[d,d]*v*v
                         for (a,b,c,d),v in riemann.items() if v!=0))
    return scalar,ric2,riem2


def verify():
    checks={}
    def zero(name,value):
        assert sp.simplify(value)==0,name
        checks[name]='exact_zero'
    t,r,u,ph=sp.symbols('t rho u phi',real=True)
    a,b,q=sp.symbols('a b q',positive=True)
    f=1+r*r;s=1-u*u
    metric=sp.diag(-a*a*f,a*a/f,b*b/s,b*b*s)
    R,Ric2,Riem2=curvature(metric,[t,r,u,ph])
    x,y=1/a**2,1/b**2
    zero('product_scalar',R-2*(y-x))
    zero('product_Ricci_squared',Ric2-2*(x*x+y*y))
    zero('product_Riemann_squared',Riem2-4*(x*x+y*y))
    Euler=sp.factor(Riem2-4*Ric2+R**2)
    Weyl2=sp.factor(Riem2-2*Ric2+R**2/3)
    zero('product_Euler',Euler+8*x*y)
    zero('product_Weyl_squared',Weyl2-sp.Rational(4,3)*(y-x)**2)
    F=sp.zeros(4);F[2,3]=q/2;F[3,2]=-q/2
    inv=metric.inv()
    F2=sp.simplify(sum(inv[i,i]*inv[j,j]*F[i,j]**2 for i in range(4) for j in range(4)))
    zero('magnetic_invariant',F2-q*q*y*y/2)

    # Independently evaluate the spin traces in the universal local a4.
    sx=sp.Matrix([[0,1],[1,0]]);sy=sp.Matrix([[0,-sp.I],[sp.I,0]]);sz=sp.diag(1,-1)
    gamma=[sp.kronecker_product(sx,sp.eye(2)),sp.kronecker_product(sy,sp.eye(2)),
           sp.kronecker_product(sz,sx),sp.kronecker_product(sz,sy)]
    XX,YY,B=sp.symbols('x y B',real=True)
    scal=2*(YY-XX);ric2=2*(XX**2+YY**2);riem2=4*(XX**2+YY**2)
    Forth=sp.zeros(4);Forth[2,3]=B;Forth[3,2]=-B
    end=-scal*sp.eye(4)/4-sp.I*sum((gamma[i]*gamma[j]*Forth[i,j]/2
                              for i in range(4) for j in range(4)),sp.zeros(4))
    omega2=0
    def Rorth(i,j,k,l):
        if all(v<2 for v in [i,j,k,l]):curv=-XX
        elif all(v>=2 for v in [i,j,k,l]):curv=YY
        else:return sp.S.Zero
        return curv*(int(i==k)*int(j==l)-int(i==l)*int(j==k))
    for i in range(4):
        for j in range(4):
            om=sum((Rorth(i,j,k,l)*gamma[k]*gamma[l]/4 for k in range(4) for l in range(4)),sp.zeros(4))-sp.I*Forth[i,j]*sp.eye(4)
            omega2+=sp.trace(om*om)
    zero('endomorphism_trace',sp.trace(end)+scal)
    zero('endomorphism_square_trace',sp.trace(end*end)-scal**2/4-4*B**2)
    zero('connection_curvature_square_trace',omega2+riem2/2+8*B**2)
    traced_a4=(60*scal*sp.trace(end)+180*sp.trace(end*end)+30*omega2
               +4*(5*scal**2-2*ric2+2*riem2))/(360*(4*sp.pi)**2)
    # Fermion determinant sign, with mostly-plus Lorentzian stress convention.
    anomaly=sp.factor(-traced_a4)
    published=(-sp.Rational(11,360)*(riem2-4*ric2+scal**2)
               +sp.Rational(1,20)*(riem2-2*ric2+scal**2/3))/(16*sp.pi**2)-2*B**2/(24*sp.pi**2)
    zero('independent_heat_kernel_vs_anomaly_coefficients',anomaly-published)
    product=sp.factor(anomaly.subs({XX:x,YY:y,B:q*y/2}))
    expected=(3*x*x+5*x*y+(3-15*q*q)*y*y)/(720*sp.pi**2)
    zero('unequal_radius_trace',product-expected)
    equal=sp.factor(product.subs(a,b))
    zero('equal_radius_trace',equal-(11-15*q*q)/(720*sp.pi**2*b**4))
    # Integer q>=1 makes 15 q^2-11 strictly positive.
    v=sp.Symbol('v',nonnegative=True)
    assert all(c>0 for c in sp.Poly(15*(v+1)**2-11,v).all_coeffs())
    checks['integer_flux_nonzero_trace']='positive_shifted_polynomial'

    # Metric calculation of the generic static spherical scalar curvature.
    SS,CC=sp.Function('S')(r),sp.Function('C')(r)
    generic=sp.diag(-SS,1/SS,CC/s,CC*s)
    generic_scalar,_,_=curvature(generic,[t,r,u,ph])
    formula=-sp.diff(SS,r,2)+2/CC-2*SS*sp.diff(CC,r,2)/CC-2*sp.diff(SS,r)*sp.diff(CC,r)/CC+SS*sp.diff(CC,r)**2/(2*CC**2)
    zero('generic_spherical_scalar',generic_scalar-formula)
    alpha=sp.Symbol('alpha',real=True)
    varphi=r*sp.atan(r)+sp.Rational(4,3)
    correction=-sp.Rational(4,3)*r*r-r*(r*r+3)*sp.atan(r)+sp.log(f)
    substituted=formula.subs({SS:f+alpha*correction,CC:1+alpha*varphi},simultaneous=True).doit()
    zero('W22_background_scalar_at_zero',substituted.subs(alpha,0))
    zero('W22_background_scalar_first_order',sp.diff(substituted,alpha).subs(alpha,0))
    T=sp.diag(r*r-2,-(r*r+4)/f**2,1/s,s)*q/(96*sp.pi**2*b*b)
    zero('W22_lifted_stress_four_trace',sp.trace(metric.subs(a,b).inv()*T))

    # The necessary trace relation fixes a ratio, not the full two pressures.
    kappa,gc=sp.symbols('kappa g_c',positive=True)
    delta=sp.Symbol('delta',real=True)
    trace_equation=2*(x-y)-kappa*expected
    leading_delta=kappa*(15*q*q-11)/(1440*sp.pi**2*b*b)
    # Loop bookkeeping: delta and kappa*T are first order; hold b fixed.
    loop=sp.Symbol('loop',real=True)
    expanded=trace_equation.subs({a:b*sp.sqrt(1+loop*leading_delta),kappa:loop*kappa},simultaneous=True)
    zero('necessary_radius_ratio_first_order',sp.diff(expanded,loop).subs(loop,0))
    re2=kappa*q*q/(8*gc*gc)
    normalized=sp.factor(leading_delta.subs(b*b,re2))
    reference_alpha=gc*gc/(4*sp.pi**2*q)
    zero('radius_shift_over_Casimir_order',normalized/reference_alpha-q/3+sp.Rational(11,45)/q)
    s0=sp.Symbol('s0',positive=True)
    zero('fixed_screening_radius_shift',normalized.subs(gc*gc,4*sp.pi**2*s0/q)
         -s0/(3*q)*(1-sp.Rational(11,15)/q**2))

    # RG-improved magnetic Maxwell density retains the anomaly under Weyl variation.
    sigma,inv_g,mu=sp.symbols('sigma inverse_g_squared mu',real=True)
    # d(1/g^2)/d log(B)=-1/(12 pi^2), so Weyl scaling B->e^-2sigma B
    # changes inverse_g by sigma/(6 pi^2).
    density=-2*B**2*(inv_g+sigma/(6*sp.pi**2))/4
    zero('running_coupling_does_not_remove_gauge_trace',sp.diff(density,sigma)+B**2/(12*sp.pi**2))
    samples=[]
    for charge in [3,8,64,256]:
        ratio=sp.Rational(charge,3)-sp.Rational(11,45*charge)
        samples.append({'q':charge,'required_radius_shift_over_alpha_exact':str(ratio),
                        'ratio_float':float(ratio)})
    return {'identities':checks,'equal_radius_anomaly':str(equal),
            'unequal_radius_anomaly':str(product),
            'required_radius_ratio_shift_first_order':str(normalized),
            'order_comparison_samples':samples,
            'inputs':'universal local heat-kernel a4; one free unit-charge Dirac determinant; no photon/graviton loops',
            'not_established':['complete renormalized stress tensor','corrected wormhole existence',
                               'physical stability or instability','novelty','global shortcut']}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_one_Dirac_loop_throat_trace_matching',
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'sympy_version':sp.__version__,'checks':verify(),'lean_certified':False,
          'necessary_background_condition_only':True,'physical_wormhole_obstruction_proved':False}
    args.output.write_bytes((json.dumps(data,indent=2)+'\n').encode())


if __name__=='__main__':main()
