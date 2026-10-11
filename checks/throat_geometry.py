"""Derive the corrected throat tensors, conditional on the drift closure.

Coordinates (tau,rho,u=cos(theta),phi), exp(p tau), r_e=1 and q=4.
Background and perturbation tensors are first derived for arbitrary S(rho)
and C(rho)=R^2, then differentiated in the MMP background correction.
This does not derive the quantum constitutive relation from the Dirac field.
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def published_radial(r, ell, screen, p, h0, h1, w0, w1, charge):
    """Appendices B and C, e=1; zeroth and first background coefficients."""
    f, a, log = 1+r*r, sp.atan(r), sp.log(1+r*r)
    dr = lambda v, n=1: sp.diff(v, r, n)
    electric = p*w1-dr(w0)
    common = 4*r*r+3*r*(f+2)*a-3*log
    base = [
        (ell+2)*h0-f*dr(h0,2)+p*f*dr(h1)+2*w0,
        -p*dr(h0)+(p*p+f*(ell+2))*h1+2*f*w1,
        p*h0-2*r*f*h1-f*f*dr(h1),
        2*ell*h0+ell*w0-f*dr(w0,2)+p*f*dr(w1)-screen*f*dr(charge),
        -2*ell*f*h1+p*dr(w0)-(p*p+ell*f)*w1+screen*p*charge,
        2*p*h0-4*r*f*h1-2*f*f*dr(h1)+p*w0-2*r*f*w1-f*f*dr(w1),
        -f*electric-p*p*charge+2*r*f*dr(charge)+f*f*dr(charge,2),
    ]
    correction = [
        -(4*ell+20+3*r*(ell+6)*a)*h0/3+common*dr(h0,2)/3
        +p*(r+f*a)*h1-p*common*dr(h1)/3-(2*f+1+2*r*f*a)*w0/f,
        p*(r+f*a)*h0/f-p*common*dr(h0)/(3*f)
        +(4*r*r*p*p-2*f*((2*ell+10)*f-3)
          +3*r*((r*r+3)*p*p-f*f*(ell+6))*a-3*p*p*log)*h1/(3*f)
        -(2*f-1+2*r*f*a)*w1,
        p*common*h0/(3*f)+f*(11*r+9*f*a)*h1/3+f*common*dr(h1)/3,
        ell*(2*f-11+12*r*a-6*log)*h0/(3*f)
        +ell*(10*f-11+6*r*(f+2)*a-6*log)*w0/(6*f)
        -(r+f*a)*dr(w0)-f*(4+3*r*a)*dr(w0,2)/3
        +p*(r+f*a)*w1+p*f*(4+3*r*a)*dr(w1)/3,
        ell*(14*f-11+12*r*(f+1)*a-6*log)*h1/3
        +p*(4+3*r*a)*dr(w0)/3
        -(8*p*p-ell*(6*f-11)+6*r*(p*p-ell*(f+2))*a+6*ell*log)*w1/6,
        p*(2*f-11+12*r*a-6*log)*h0/(3*f)
        +4*f*(10*r+(9*f-3)*a)*h1/3
        +f*(14*f-11+12*r*(f+1)*a-6*log)*dr(h1)/3
        +p*(10*f-11+6*r*(f+2)*a-6*log)*w0/(6*f)
        +3*f*(r+f*a)*w1+f*(6*f-11+6*r*(f+2)*a-6*log)*dr(w1)/6,
        -p*p*common*charge/(3*f)-f*(11*r+9*f*a)*dr(charge)/3
        -f*common*dr(charge,2)/3,
    ]
    return base, correction


def derive_geometry():
    r, u, ell, screen, p = sp.symbols('rho u L s p', real=True)
    h0, h1, w0, w1, charge = [sp.Function(n)(r) for n in ['h0','h1','w0','w1','X']]
    sf, cf = sp.Function('S')(r), sp.Function('C')(r)
    y = sp.Function('Y')(u)
    yp = sp.diff(y,u)
    y2 = (2*u*yp-ell*y)/(1-u*u)
    angular = {sp.diff(y,u,3):sp.diff(y2,u).subs(sp.diff(y,u,2),y2),sp.diff(y,u,2):y2}

    def simplify(v):
        return sp.factor(v.xreplace(angular).subs(sp.diff(y,u,2),y2))

    def derivative(v, index, pert=False):
        return (p*v if pert else sp.S.Zero) if index==0 else (sp.diff(v,r) if index==1 else (sp.diff(v,u) if index==2 else sp.S.Zero))

    metric = sp.diag(-sf,1/sf,cf/(1-u*u),cf*(1-u*u))
    inverse = sp.diag(-1/sf,sf,(1-u*u)/cf,1/(cf*(1-u*u)))
    h = sp.zeros(4)
    h[0,3]=h[3,0]=-h0*(1-u*u)*yp
    h[1,3]=h[3,1]=-h1*(1-u*u)*yp
    dinverse=-inverse*h*inverse
    assert sp.trace(inverse*h)==0
    connection, dconnection = {}, {}
    for a in range(4):
        for b in range(4):
            for c in range(b,4):
                v=dv=sp.S.Zero
                for d in range(4):
                    dg=derivative(metric[d,c],b)+derivative(metric[d,b],c)-derivative(metric[b,c],d)
                    dh=derivative(h[d,c],b,True)+derivative(h[d,b],c,True)-derivative(h[b,c],d,True)
                    v+=inverse[a,d]*dg/2
                    dv+=(inverse[a,d]*dh+dinverse[a,d]*dg)/2
                connection[a,b,c]=connection[a,c,b]=sp.factor(v)
                dconnection[a,b,c]=dconnection[a,c,b]=simplify(dv)
    ricci,dricci=sp.zeros(4),sp.zeros(4)
    for a in range(4):
        for b in range(a,4):
            v=dv=sp.S.Zero
            for c in range(4):
                v+=derivative(connection[c,a,b],c)-derivative(connection[c,a,c],b)
                dv+=derivative(dconnection[c,a,b],c,True)-derivative(dconnection[c,a,c],b,True)
                for d in range(4):
                    v+=connection[c,d,c]*connection[d,a,b]-connection[c,d,b]*connection[d,a,c]
                    dv+=dconnection[c,d,c]*connection[d,a,b]+connection[c,d,c]*dconnection[d,a,b]
                    dv-=dconnection[c,d,b]*connection[d,a,c]+connection[c,d,b]*dconnection[d,a,c]
            ricci[a,b]=ricci[b,a]=sp.factor(v)
            dricci[a,b]=dricci[b,a]=simplify(dv)
    scalar=sp.factor(sp.trace(inverse*ricci))
    assert simplify(sp.trace(inverse*dricci)+sp.trace(dinverse*ricci))==0
    field=sp.zeros(4);field[2,3]=2;field[3,2]=-2
    electric=p*w1-sp.diff(w0,r)
    dfield=sp.zeros(4)
    dfield[0,1]=electric*y;dfield[0,2]=-w0*yp;dfield[1,2]=-w1*yp
    dfield=dfield-dfield.T
    field_up=inverse*field*inverse
    dfield_up=dinverse*field*inverse+inverse*dfield*inverse+inverse*field*dinverse
    invariant=sum(field[a,b]*field_up[a,b] for a in range(4) for b in range(4))
    dinvariant=sum(dfield[a,b]*field_up[a,b]+field[a,b]*dfield_up[a,b] for a in range(4) for b in range(4))
    stress=field*inverse*field.T-metric*invariant/4
    dstress=dfield*inverse*field.T+field*dinverse*field.T+field*inverse*dfield.T-h*invariant/4-metric*dinvariant/4
    base_einstein=ricci-metric*scalar/2-stress/2
    einstein=(dricci-h*scalar/2-dstress/2).applyfunc(simplify)
    maxwell=sp.Matrix([simplify(sum(derivative(cf*dfield_up[a,b],a,True) for a in range(4))) for b in range(4)])
    print('general_spherical_tensors_derived',flush=True)
    return ((r,u,ell,screen,p), (h0,h1,w0,w1,charge), (sf,cf), (y,yp),
            (simplify,derivative), (metric,inverse,h,dinverse),
            (connection,dconnection), (field,field_up),
            (base_einstein,einstein,maxwell))


def verify():
    ((r,u,ell,screen,p), (h0,h1,w0,w1,charge), (sf,cf), (y,yp),
     (simplify,derivative), (metric,inverse,h,dinverse),
     (connection,dconnection), (field,field_up),
     (base_einstein,einstein,maxwell))=derive_geometry()
    electric=p*w1-sp.diff(w0,r)

    f=1+r*r
    gamma=-4*r*r/3-r*(r*r+3)*sp.atan(r)+sp.log(f)
    phi=r*sp.atan(r)+sp.Rational(4,3)
    jet_rules, base_jets, tangent_jets={},{},{}
    for func,base,corr in [(sf,f,gamma),(cf,sp.S.One,phi)]:
        for n in range(3):
            j=sp.Symbol(str(func.func)+'_'+str(n))
            jet_rules[sp.diff(func,r,n)]=j
            base_jets[j]=sp.diff(base,r,n)
            tangent_jets[j]=sp.diff(corr,r,n)

    def expansion(v):
        polynomial=v.xreplace(jet_rules)
        base=polynomial.xreplace(base_jets)
        corr=sum(sp.diff(polynomial,j).xreplace(base_jets)*d for j,d in tangent_jets.items())
        return simplify(base),simplify(corr)

    def zeroth(v):return expansion(v)[0]
    g0=metric.applyfunc(zeroth);inv0=inverse.applyfunc(zeroth);dinv0=dinverse.applyfunc(zeroth)
    tq=sp.diag(r*r-2,-(r*r+4)/f**2,1/(1-u*u),1-u*u)/3
    dtq=sp.zeros(4)
    kval=[(1-2*r*r)/f,-(2*r*r+5)/f]
    for a,ha,wa in [(0,h0,w0),(1,h1,w1)]:
        tilt=-(ha+wa/2)*(1-u*u)*yp
        dtq[a,3]=dtq[3,a]=(h[a,3]+kval[a]*tilt)/3
    coeff=[inv0[0,0]*kval[0]*(h0+w0/2),inv0[1,1]*kval[1]*(h1+w1/2)]
    jpol=sp.Matrix([-ell*coeff[0]*y/3,-ell*coeff[1]*y/3,
                    -(1-u*u)*(p*coeff[0]+sp.diff(coeff[1],r))*yp/3,0])
    jsch=sp.Matrix([screen*sp.diff(charge,r)*y,-screen*p*charge*y,0,0])
    assert simplify(sum(derivative(jpol[a],a,True) for a in range(4)))==0
    assert simplify(sum(derivative(jsch[a],a,True) for a in range(4)))==0

    checks={}
    def zero(name,v):
        result=simplify(v)
        assert result==0,(name,result)
        checks[name]='exact_zero'

    for a in range(4):
        for b in range(a,4):
            b0,b1=expansion(base_einstein[a,b])
            zero('background_Einstein_0_'+str(a)+str(b),b0)
            zero('background_Einstein_1_'+str(a)+str(b),b1-tq[a,b])
    for a in range(4):
        zero('background_Maxwell_'+str(a),sum(derivative(cf*field_up[b,a],b) for b in range(4)))
    mixed=inv0*tq
    dmixed=inv0*dtq+dinv0*tq
    for nu in range(4):
        div=sum(derivative(dmixed[mu,nu],mu,True) for mu in range(4))
        for mu in range(4):
            for lam in range(4):
                div+=zeroth(connection[mu,mu,lam])*dmixed[lam,nu]-zeroth(connection[lam,mu,nu])*dmixed[mu,lam]
                div+=zeroth(dconnection[mu,mu,lam])*mixed[lam,nu]-zeroth(dconnection[lam,mu,nu])*mixed[mu,lam]
        zero('quantum_force_balance_'+str(nu),div+sum(field[nu,lam]*jpol[lam]/2 for lam in range(4)))
    print('background_and_quantum_Ward_checks_passed',flush=True)

    e0,e1={},{}
    for a in range(4):
        for b in range(a,4):
            e0[a,b],v=expansion(einstein[a,b]);e1[a,b]=simplify(v-dtq[a,b])
            if (a,b) not in [(0,3),(1,3),(2,3)]:
                zero('other_Einstein_0_'+str(a)+str(b),e0[a,b]);zero('other_Einstein_1_'+str(a)+str(b),e1[a,b])
    m0,m1=[],[]
    for b in range(4):
        a,c=expansion(maxwell[b]);m0.append(simplify(a-jsch[b]));m1.append(simplify(c-jpol[b]))
    zero('Maxwell_phi_0',m0[3]);zero('Maxwell_phi_1',m1[3])
    xeq=-electric-p*p*charge/sf+sp.diff(sf*sp.diff(charge,r),r)
    x0,x1=expansion(xeq)
    actual0=[e0[0,3],e0[1,3],e0[2,3],m0[0],m0[1],m0[2],x0]
    actual1=[e1[0,3],e1[1,3],e1[2,3],m1[0],m1[1],m1[2],x1]
    names=['Einstein_tau_phi','Einstein_rho_phi','Einstein_angular','Maxwell_tau','Maxwell_rho','Maxwell_angular','charge_anomaly']
    paper0,paper1=published_radial(r,ell,screen,p,h0,h1,w0,w1,charge)
    ratios=[-(1-u*u)*yp/2,-(1-u*u)*yp/(2*f),(-ell*y+2*u*yp)/(2*f),y/f,y,(1-u*u)*yp/f,1/f]
    for i,name in enumerate(names):
        zero(name+'_order_0',actual0[i]-ratios[i]*paper0[i])
        zero(name+'_order_1',actual1[i]-ratios[i]*paper1[i])
        print(name,'both_orders_exact_zero',flush=True)
    return {'identities':checks,'charge_conservation':'exact_zero_for_Schwinger_and_polarisation_currents',
            'radial_matches':dict.fromkeys(names,'both_orders_exact_zero'),
            'independent_background_and_linear_tensor_derivation':True,
            'background_order':1,'perturbation_order':1,'quantum_drift_closure_assumed':True}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_independent_throat_tensor_derivation',
          'source':'https://arxiv.org/abs/2610.09847v1',
          'source_equations':['2.5-2.6','5.7-5.8','6.6','6.21-6.24','B.1-B.3','B.6','B.12-B.14','C.1-C.8'],
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'sympy_version':sp.__version__,'checks':verify(),
          'quantum_closure_derived_from_Dirac_theory':False,'global_matching_verified':False,
          'lean_certified':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
