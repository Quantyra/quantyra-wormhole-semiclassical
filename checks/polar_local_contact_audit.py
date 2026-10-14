"""S074 audit of a preserved contact draft and its missing closure obligations.

Algebra through finite_checks is adapted from the preexisting untracked
polar_dipole_contacts.py, SHA-256
14ca4c38524b91cf231ef50d50df29f59277cca99614241dd2289196c278427a.
The original is preserved. No full physical determinant is asserted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
import scipy
from scipy.integrate import quad
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def zero(x):
    reduced=s.simplify(s.expand(x))
    assert reduced==0,str(reduced)


def curvature_jet():
    # Four-dimensional metric jets at w=0, partial_x w=k. No KK formula used.
    A,B,H,R,k=s.symbols('A B H R k',positive=True)
    ap,bp,hp,rp,app,bpp,hpp,rpp=s.symbols('ap bp hp rp app bpp hpp rpp')
    g=s.diag(-A,B,H,R);gi=g.inv()
    dg=[s.zeros(4) for _ in range(4)]
    dg[1][0,3]=dg[1][3,0]=-R*k
    dg[2]=s.diag(-ap,bp,hp,rp)
    ddg=[[s.zeros(4) for _ in range(4)] for _ in range(4)]
    ddg[1][1][0,0]=2*R*k*k
    ddg[1][2][0,3]=ddg[1][2][3,0]=-rp*k
    ddg[2][1]=ddg[1][2].copy();ddg[2][2]=s.diag(-app,bpp,hpp,rpp)
    dgi=[-gi*x*gi for x in dg]
    gam={};dgam={}
    for i in range(4):
        for j in range(4):
            for q in range(4):
                gam[i,j,q]=sum(gi[i,z]*(dg[j][z,q]+dg[q][z,j]-dg[z][j,q])/2 for z in range(4))
                for v in range(4):
                    dgam[v,i,j,q]=sum((dgi[v][i,z]*(dg[j][z,q]+dg[q][z,j]-dg[z][j,q])+gi[i,z]*(ddg[v][j][z,q]+ddg[v][q][z,j]-ddg[v][z][j,q]))/2 for z in range(4))
    ric=s.zeros(4)
    for i in range(4):
        for j in range(4):
            ric[i,j]=sum(dgam[q,q,i,j]-dgam[j,q,i,q]+sum(gam[q,i,j]*gam[z,q,z]-gam[z,i,q]*gam[q,j,z] for z in range(4)) for q in range(4))
    scalar=s.expand(s.trace(gi*ric));coefficient=scalar.coeff(k,2)
    zero(coefficient-R/(2*A*B))
    # Original Maxwell contraction at the same horizontal-frame point.
    F,G,c,aa,V,ht,hx=s.symbols('F G c aa V ht hx')
    em=s.zeros(4);em[0,1]=c*G+aa*V*F;em[1,0]=-em[0,1]
    em[2,0]=ht;em[0,2]=-ht;em[2,1]=hx;em[1,2]=-hx
    fsq=sum(em[i,j]*em[u,v]*gi[i,u]*gi[j,v] for i in range(4) for j in range(4) for u in range(4) for v in range(4))
    zero(fsq-(-2*(c*G+aa*V*F)**2/(A*B)+2/H*(-ht*ht/A+hx*hx/B)))
    return {'Ricci_F_squared_coefficient':str(coefficient),'Maxwell_horizontal_and_mixed_contractions':'exact','scope':'Arbitrary diagonal axisymmetric background metric and its first/second c derivatives; connection-curvature term independently reconstructed.'}


def contact_algebra():
    c,chi,a,J,kg,ke,P,mu,t=s.symbols('c chi a J kg ke P mu t',positive=True)
    Y=s.legendre(2,c);V=(1-c*c)*s.diff(Y,c)
    integ=lambda f:s.integrate(s.expand(f),(c,-1,1))
    grav=3*s.Rational(1,4)*integ((1-c*c)*(1+4*t*chi*Y+t*t*(8*chi*chi+J/2)*Y*Y))
    electric=s.Rational(3,2)*integ(c*c*(1+2*t*chi*Y+t*t*(2*chi*chi+J/2)*Y*Y))
    cross=s.Rational(3,2)*integ(c*V*t*a*(1+2*t*chi*Y))
    extra=s.Rational(3,2)*integ(V*V*t*t*a*a)
    zero(grav-(1-4*t*chi/5+t*t*(8*chi*chi/7+J/14)))
    zero(electric-(1+4*t*chi/5+t*t*(22*chi*chi/35+11*J/70)))
    zero(cross-(6*t*a/5+12*t*t*chi*a/35));zero(extra-72*t*t*a*a/35)
    D=s.diag(kg,ke)
    K1=s.Matrix([[-4*kg*chi/5,6*ke*P*a/5],[6*ke*P*a/5,4*ke*chi/5]])
    K2=s.Matrix([[kg*(8*chi*chi/7+J/14)+ke*P*P*72*a*a/35,ke*P*12*chi*a/35],[ke*P*12*chi*a/35,ke*(22*chi*chi/35+11*J/70)]])
    inv=D.inv();inverse=inv-t*inv*K1*inv+t*t*(inv*K1*inv*K1*inv-inv*K2*inv)
    v=s.Matrix([-P,1]);mass=s.expand(mu*(v.T*inverse*v)[0])
    r=mu*P*P/kg;d=mu/ke
    M1=4*(r-d)*chi/5+12*r*a/5
    M2=(-88*r+2*d)*chi*chi/175+24*r*chi*a/35+r*(s.Rational(36,25)-s.Rational(108,175)*r/d)*a*a-(5*r+11*d)*J/70
    zero(mass-(r+d+t*M1+t*t*M2))
    # Lorentzian tensor inversion, including an off-diagonal base perturbation.
    p,q=s.symbols('p q');eta=s.diag(-1,1);h=s.Matrix([[p,q],[q,p]])
    hc=eta*h*eta;j=p*p-q*q;eps=s.Matrix([[0,-1],[1,0]])
    U1=hc/5;U2=j*eta/14
    Ui=eta-t*eta*U1*eta+t*t*(eta*U1*eta*U1*eta-eta*U2*eta)
    Z=-eps*Ui*eps.T
    assert (Z-(eta+t*hc/5-11*t*t*j*eta/350)).applyfunc(s.factor)==s.zeros(2)
    f=s.series(s.sqrt((Z*eta).det()),t,0,3).removeO()
    zero(f-(1-9*t*t*j/175))
    heatM2=s.expand(M2+9*(r+d)*J/175)
    expected=M2+(5*r+11*d)*J/70-(7*r+37*d)*J/350
    zero(heatM2-expected)
    zero(c*Y-s.Rational(2,5)*s.legendre(1,c)-s.Rational(3,5)*s.legendre(3,c))
    assert integ(s.legendre(3,c)*c*Y)==s.Rational(6,35)
    return {'K_gravity_series':str(grav),'K_electric_series':str(electric),'K_cross_series':str(cross),'K_extra_series':str(extra),'mass_linear':str(s.factor(M1)),'mass_quadratic':str(s.factor(M2)),'kinetic_seagull':'-11 J g^ab/350','sqrt_det_Z':'1-9 J/175','bulk_log_mass_quadratic':str(s.factor(heatM2)),'omitted_l3_overlap':'6/35'}


def cubic_matching():
    C,b,G,e=s.symbols('C b G e',positive=True);kap=8*s.pi*G
    W=s.Matrix([[2*b*b/kap*(1+C/3),2*b*b*e],[2*b*b*e,6]])
    R=s.Matrix([3*G*(C-1)/b**4,9*G*e/b**2])
    back=lambda x:s.factor(s.expand(x).subs(e*e,(C+1)/(kap*b*b)))
    direct=s.Matrix([-2*(C*C+11*C+6)/(15*b*b),-16*s.pi*G*e*(C+8)/(5*b*b)])
    assert (direct+s.Rational(8,15)*s.pi*W*R).applyfunc(back)==s.zeros(2,1)
    # Divergence identity behind the integration by parts: free scalar stress.
    # Check it independently in a curved conformal metric using arbitrary jets.
    x=s.symbols('x',real=True);w,m2=s.symbols('w m2',real=True)
    B=1+x*x;f=s.Function('f')(x)
    # phi=e^-iwt f(x), (Box-m2)phi=0.
    fp=s.diff(f,x);fpp=(m2-w*w/B)*f/B-s.diff(B,x)*fp/B
    # D_a[(D^a phi D^b phi)-g^ab(Dphi)^2/2] = m2 phi D^b phi.
    grad2=w*w*f*f/B+B*fp*fp
    stressTT=-w*w*f*f/B**2+grad2/(2*B)
    stressTx=s.I*w*f*fp
    stressxx=B*B*fp*fp-B*grad2/2
    divT=-2*s.I*w*stressTT+s.diff(stressTx,x)+s.diff(B,x)/B*stressTx
    divx=-2*s.I*w*stressTx+s.diff(stressxx,x)+B*s.diff(B,x)*stressTT/2-s.diff(B,x)*stressxx/(2*B)
    zero(s.expand(divT-s.I*w*m2*f*f/B).subs(s.diff(f,x,2),fpp))
    zero(s.expand(divx-m2*B*f*fp).subs(s.diff(f,x,2),fpp))
    return {'on_shell_coefficients':[str(x) for x in direct],'harmonic_conversion':'8pi/15','curved_stress_divergence':'exact for arbitrary radial data satisfying the free wave equation'}


def exact_operator(C,chi,a,p):
    Y=lambda c:(3*c*c-1)/2
    V=lambda c:3*c*(1-c*c)
    d=lambda c:math.sqrt(1-p*p*Y(c)**2)
    I=lambda f:quad(f,-1,1,epsabs=2e-13,epsrel=2e-13)[0]
    ke=C+1;mu=2*(C+1)/C
    K=np.array([[.75*I(lambda c:(1-c*c)*math.exp(4*chi*Y(c))/d(c))+1.5*ke*a*a*I(lambda c:V(c)**2*math.exp(2*chi*Y(c))/d(c)),1.5*ke*a*I(lambda c:c*V(c)*math.exp(2*chi*Y(c))/d(c))],[0,1.5*ke*I(lambda c:c*c*math.exp(2*chi*Y(c))/d(c))]])
    K[1,0]=K[0,1]
    assert np.linalg.eigvalsh(K)[0]>0
    v=np.array([-1.,1.]);M=mu*float(v@np.linalg.solve(K,v))
    U=[.75*I(lambda c:(1-c*c)*(1-sign*p*Y(c))/d(c)) for sign in [1,-1]]
    Z=np.array([1/U[1],1/U[0]])
    return M,Z


def finite_checks():
    rows=[];heat=[]
    for C in [1.,4.,25.]:
        r=2*(C+1)/C;d=2/C
        for shape in [(0.7,0.4,0.9),(-0.6,0.2,0.5)]:
            errors=[]
            for t in [.02,.01,.005]:
                chi,a,p=np.array(shape)*t;J=p*p
                M,Z=exact_operator(C,chi,a,p)
                Mpoly=r+d+4*(r-d)*chi/5+12*r*a/5+(-88*r+2*d)*chi*chi/175+24*r*chi*a/35+r*(36/25-108/175*r/d)*a*a-(5*r+11*d)*J/70
                Zpoly=np.array([1+p/5-11*J/350,1-p/5-11*J/350])
                error=max(abs(M-Mpoly)/(r+d),float(np.max(abs(Z-Zpoly))))
                errors.append(error)
            assert errors[0]/errors[1]>5 and errors[1]/errors[2]>5
            rows.append({'C':C,'shape':list(shape),'amplitudes':[.02,.01,.005],'relative_operator_remainder':errors,'halving_ratios':[errors[0]/errors[1],errors[1]/errors[2]]})
        M,Z=exact_operator(C,.014,.008,.018);f=math.sqrt(float(np.prod(Z)))
        n=np.arange(-220,221,dtype=float)
        coeff=[];maxerr=0.
        for t in [.02,.01]:
            density=math.exp(-t*M)*float(np.sum(np.exp(-t*Z[0]*n*n))*np.sum(np.exp(-t*Z[1]*n*n)))/(4*math.pi**2)
            expected=math.exp(-t*M)/(4*math.pi*t*f)
            rel=abs(density/expected-1);assert rel<2e-12;maxerr=max(maxerr,rel)
            coeff.append(4*math.pi*density-1/(t*f))
        extrap=2*coeff[1]-coeff[0];bound=.02**2*M**3/(4*f)
        assert abs(extrap+M/f)<bound+1e-10
        heat.append({'C':C,'Fourier_cutoff':220,'max_heat_trace_relative_error':maxerr,'extrapolated_heat_coefficient':extrap,'analytic_heat_coefficient':-M/f,'explicit_extrapolation_error_bound':bound})
    return {'operator_comparisons':rows,'heat_trace_checks':heat}


def closure_checks():
    c,t,aa,bb,v=s.symbols('c t aa bb v',real=True)
    P1=s.legendre(1,c);P2=s.legendre(2,c);P3=s.legendre(3,c)
    zero(P1*P2-s.Rational(2,5)*P1-s.Rational(3,5)*P3)
    overlap=s.integrate(P1*P2*P3,(c,-1,1))
    zero(overlap-s.Rational(6,35))
    normalized=s.sqrt(21)*overlap/2
    zero(normalized**2-s.Rational(27,175))
    total_norm=s.Rational(3,2)*s.integrate((P1*P2)**2,(c,-1,1))
    zero(total_norm-s.Rational(11,35))
    # Exact finite block determinant, independent of any continuum truncation.
    H=s.Matrix([[aa,t*v],[t*v,bb]])
    det=H.det()
    zero(det-(aa*bb-t*t*v*v))
    log2=s.diff(s.log(det),t,2).subs(t,0)/4
    zero(log2+v*v/(2*aa*bb))
    # Eliminating a genuinely invertible constraint adds its quartic Schur term.
    q,J,A=s.symbols('q J A',real=True)
    qsol=-J/A
    zero((A*q*q/2+q*J).subs(q,qsol)+J*J/(2*A))
    diagnostics=[]
    for a,b,coupling in [(2.,5.,.3),(1.,12.,.8),(3.,4.,-.5)]:
        errors=[]
        for step in [.02,.01,.005]:
            sign,ld=np.linalg.slogdet([[a,step*coupling],[step*coupling,b]])
            assert sign>0
            exact=.5*(ld-math.log(a*b))
            leading=-step*step*coupling*coupling/(2*a*b)
            remainder=abs(exact-leading)
            x=step*step*coupling*coupling/(a*b)
            # |(log(1-x)+x)/2| <= x^2/[4(1-x)].
            bound=x*x/(4*(1-x))
            assert remainder<=bound+8e-16
            errors.append({'step':step,'remainder':remainder,'analytic_bound':bound})
        diagnostics.append({'a':a,'b':b,'v':coupling,'cases':errors})
    # Kinematic hierarchy, not a bound on the constrained mixed vertex.
    C=s.symbols('C',positive=True)
    mu3=120/(C*(C+12+s.sqrt(C*C+24*C+24)))
    assert s.limit(C*C*mu3,C,s.oo)==60
    hierarchy=[]
    for k in [2000,4000,10000]:
        cc=math.sqrt(1+12*math.pi*k**6)
        light3=120/(cc*(cc+12+math.sqrt(cc*cc+24*cc+24)))
        dipole=2+4/cc
        hierarchy.append({'k':k,'C':cc,'l3_light_mass_squared_a_squared':light3,
                          'dipole_mass_squared_a_squared':dipole,
                          'light_to_dipole_mass_squared_ratio':light3/dipole})
    return {'unnormalized_l1_l2_l3_overlap':str(overlap),
            'normalized_overlap_squared':str(s.factor(normalized**2)),
            'l1_only_retains_norm_squared':'4/25',
            'total_multiplication_norm_squared':str(total_norm),
            'mixed_block_logdet_quadratic':str(log2),
            'constraint_schur_quartic':'-J^2/(2 A), only on an invertible constrained domain',
            'determinant_diagnostics':diagnostics,
            'l3_light_asymptotic':'mu_l3_minus = 60/C^2 + o(C^-2)',
            'hierarchy':hierarchy,
            'scope':'Angular overlap and block algebra expose missing terms; no physical mixed-vertex coefficient, full determinant or loop-error size is supplied.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out={'story':'S074','scope':'Aligned axisymmetric dipole block on a specified finite polar background; direct cubic matching and unreduced local quartic contacts, with scalar-block bulk logarithm only.',
         'four_dimensional_tensor':curvature_jet(),'contacts':contact_algebra(),'cubic_matching':cubic_matching(),'numerical':finite_checks(),'closure':closure_checks()}
    out['versions']={'sympy':s.__version__,'numpy':np.__version__,'scipy':scipy.__version__}
    out['sources']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['checks/polar_local_contact_audit.py','notes/polar-local-contact-audit.md','notes/dipole-quadrupole-source.md','notes/polar-canonical-vacuum.md','notes/polar-dipole-dressing.md','notes/quadrupole-product-response.md']}
    out['status']='PASS';out['goal_complete']=False;out['S074_complete']=False
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print('PASS: aligned contact audit, harmonic omission and constraint/determinant closure checks; S074 remains incomplete',flush=True)


if __name__=='__main__':main()
