"""S145: exact general-jet loop coefficients and quantitative interior response.

Numerical quadratures and finite matrices are diagnostics, not continuum proof.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss
from scipy.special import kv
import scipy

ROOT=Path(__file__).resolve().parents[1]
STEM="general-interior-dirac-response"
INPUTS=[
    "checks/general_interior_dirac_response.py",
    "notes/general-interior-dirac-response.md",
    "notes/coincident-normal-dirac-response.md",
    "notes/curved-remote-dirac-response.md",
    "notes/flat-variable-dirac-three-point.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-curved-quantum-tt.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS=[]
NUMERICS=[]
MOMENTUM_CASES=[]


def check(name,passed,detail=""):
    if any(row["name"]==name for row in ROWS):
        raise ValueError("duplicate check label")
    row={"name":name,"passed":bool(passed),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def exact_loop_checks():
    x,y=s.symbols("x y",real=True)
    den=1+x*x+y*y
    gx=s.Matrix([[0,1],[1,0]])
    gy=s.Matrix([[0,-s.I],[s.I,0]])
    N=s.eye(2)-s.I*(x*gx+y*gy)
    # Each tuple is the actual coordinate monomial coefficient in h.
    mono={
        "x":(1,0,s.Integer(1)), "qx":(2,0,s.Rational(1,2)),
        "qy":(0,2,s.Rational(1,2)), "qxy":(1,1,s.Integer(1)),
        "cx":(3,0,s.Rational(1,6)), "cxyy":(1,2,s.Rational(1,2)),
        "fourx":(4,0,s.Rational(1,24)),
    }

    def differentiate(P,p,var):
        return (P.diff(var)*den-2*p*var*P).applyfunc(s.expand),p+1

    @lru_cache(maxsize=None)
    def ordered(labels):
        P,p=N,1
        for label in reversed(labels):
            dx,dy,coef=mono[label]
            for _ in range(dx):
                P,p=differentiate(P,p,x)
            for _ in range(dy):
                P,p=differentiate(P,p,y)
            P=(coef*s.I**(dx+dy)*N*P).applyfunc(s.expand)
            p+=1
        return s.expand((-1)**(len(labels)+1)*s.trace(P)),p

    def integrate_polynomial(poly,p):
        result=s.Integer(0)
        for (ix,iy),coef in s.Poly(poly,x,y).terms():
            if ix%2 or iy%2:
                continue
            i,j=ix//2,iy//2
            assert p-i-j-1>0
            result+=coef*s.gamma(i+s.Rational(1,2))*s.gamma(j+s.Rational(1,2))*\
                s.gamma(p-i-j-1)/(4*s.pi**2*s.gamma(p))
        return s.simplify(result)

    cases=[
        ("one_quadratic",("qx",),-s.Rational(1,12)/s.pi),
        ("one_quartic",("fourx",),-s.Rational(1,120)/s.pi),
        ("two_linear",("x","x"),s.Rational(1,12)/s.pi),
        ("two_quadratic_x",("qx","qx"),s.Rational(1,20)/s.pi),
        ("two_quadratic_cross",("qx","qy"),s.Rational(1,15)/s.pi),
        ("two_quadratic_shear",("qxy","qxy"),s.Rational(1,30)/s.pi),
        ("linear_cubic_x",("x","cx"),s.Rational(1,15)/s.pi),
        ("linear_cubic_xyy",("x","cxyy"),s.Rational(1,15)/s.pi),
        ("three_x_x_qx",("x","x","qx"),-s.Rational(17,60)/s.pi),
        ("three_x_x_qy",("x","x","qy"),-s.Rational(3,20)/s.pi),
        ("three_x_x_qxy",("x","x","qxy"),s.Integer(0)),
        ("four_linear",("x","x","x","x"),s.Rational(7,40)/s.pi),
        ("odd_linear_quadratic",("x","qx"),s.Integer(0)),
        ("odd_three_linear",("x","x","x"),s.Integer(0)),
    ]
    for name,labels,target in cases:
        items=[ordered(order) for order in sorted(set(permutations(labels)))]
        assert len(set(p for _,p in items))==1
        poly=s.expand(sum(poly for poly,_ in items))
        p=items[0][1]
        value=integrate_polynomial(poly,p)
        check("exact_loop_"+name,s.simplify(value-target)==0,str(value))
        if name in ("two_linear","linear_cubic_x","linear_cubic_xyy",
                    "three_x_x_qx","three_x_x_qy","four_linear"):
            MOMENTUM_CASES.append((name,s.lambdify((x,y),poly/den**p,"numpy"),float(target)))


def jet_and_action_checks():
    x,y=s.symbols("x y",real=True)
    lx,ly,a,b,c,tx,txy,tyx,ty,B=s.symbols("lx ly a b c tx txy tyx ty B",real=True)
    h=lx*x+ly*y+(a*x*x+2*c*x*y+b*y*y)/2+\
        tx*x**3/6+txy*x*y*y/2+tyx*x*x*y/2+ty*y**3/6+B*x**4/24
    # Fourth Taylor jet of log(1+h), sufficient for the exact point identities.
    phi=s.Integer(0)
    for power in range(1,5):
        expanded=s.Poly(s.expand(h**power),x,y)
        truncated=sum(coef*x**i*y**j for (i,j),coef in expanded.terms() if i+j<=4)
        phi+=s.Rational((-1)**(power+1),power)*truncated
    lap=lambda f:s.diff(f,x,2)+s.diff(f,y,2)
    at0=lambda f:s.expand(f.subs({x:0,y:0}))
    grad=lambda f:s.Matrix([s.diff(f,x),s.diff(f,y)])
    A=at0(lap(phi))
    Q=at0(lap(lap(phi))-3*lap(phi)**2-4*grad(phi).dot(grad(lap(phi)))
          +4*grad(phi).dot(grad(phi))*lap(phi))
    ll=lx*lx+ly*ly
    H=s.Matrix([[a,c],[c,b]])
    L=s.Matrix([lx,ly])
    expected=B-4*(a+b)**2-2*s.trace(H*H)-8*(lx*(tx+txy)+ly*(tyx+ty))+\
        18*(a+b)*ll+16*(L.T*H*L)[0]-21*ll**2
    check("logarithmic_second_jet",s.expand(A-(a+b-ll))==0)
    check("logarithmic_fourth_invariant",s.expand(Q-expected)==0)

    # Direct Euler-Lagrange formula, without normal-point restrictions.
    f=s.Function("f")(x,y)
    alpha=s.symbols("alpha")
    source=-8*alpha*s.exp(-2*f)*lap(f)**2+8*alpha*lap(s.exp(-2*f)*lap(f))
    q=lap(lap(f))-3*lap(f)**2-4*grad(f).dot(grad(lap(f)))+\
        4*grad(f).dot(grad(f))*lap(f)
    check("general_covariant_curvature_density",s.simplify(source-8*alpha*s.exp(-2*f)*q)==0)
    check("scalar_area_factor",s.simplify(s.exp(-2*f)*source-8*alpha*s.exp(-4*f)*q)==0)
    M,fc=s.symbols("M fc",positive=True)
    check("local_mass_restores_covariant_coefficient",
          s.simplify(1/(M*s.exp(fc))**2-s.exp(-2*fc)/M**2)==0)
    check("proper_flat_vacuum_prefactor",
          s.simplify(M*s.exp(-fc)*(M*s.exp(fc))-M*M)==0)


def rational_bound_checks():
    rm=F(102,100)
    dc=F(11,1000)
    check("reference_mass_real_part",(1-dc)*(1-dc*dc/2)>F(98,100))
    check("reference_mass_modulus",1/(1-dc)<rm)
    check("reference_mass_phase_ratio",1/(1-dc*dc/2)<rm)
    check("scalar_area_exponential",1/(1-2*dc)<F(103,100))
    check("local_gradient_coefficient",F(10,7)*dc<F(2,100))
    linear=(8/(32*(1-F(26,1000)))+F(16,1000)+dc*64)/(1-F(5,100))
    check("global_linear_potential_bound",linear<2,float(linear))
    quadratic=(8/(1-F(26,1000))+F(16,1000)*64+dc*64**2)/(1-F(5,100))+F(2,100)*64
    check("global_quadratic_remainder",quadratic<100,float(quadratic))
    check("uniform_small_potential",F(5,100)/(1-F(5,100))<F(6,100))
    check("global_cubic_remainder",F(202,100)*64**2+9*64<2*10**4)
    check("global_fourth_remainder",F(202,100)*64**3+9*64**2+2*10**4*64<10**8)
    check("global_fifth_remainder",F(202,100)*64**4+9*64**3+2*10**4*64**2+10**6*64<10**12)
    cmin=F(98,100)*F(7,10)-F(1,10)
    check("complex_kernel_L1",2*F(10,7)*rm/cmin<5)
    check("complex_kernel_L_three_halves",8*F(10,7)*rm**3/(F(3,2)*cmin)<14)
    check("complex_two_edge_fifth",16*rm**2/F(98,100)**5<20)
    p={1:F(2,100),2:F(9),3:F(2*10**4),4:F(10**6)}
    B={1:F(2),2:F(100),3:F(2*10**4),4:F(10**8),5:F(10**12)}

    @lru_cache(None)
    def remainder(N,d):
        if N==1:
            return B[d+1]
        if d==N-1:
            return F(2)**N
        return B[d-N+2]*2**(N-1)+sum(p[a]*remainder(N-1,d-a)
                                     for a in range(1,d-N+2))

    targets=[F(10**12),F(204180000),F(82644),F(8080808,10000)]
    coeff=[remainder(n,4) for n in range(1,5)]
    check("taylor_product_recursion",coeff==targets,[str(v) for v in coeff])
    pref=28*math.factorial(5)*10**5
    one=20*rm**2*coeff[0]
    low=sum(pref*5**(n-2)*rm**(n+1)*coeff[n-1] for n in range(2,5))
    high=pref*2**5*5**3*rm**6/(1-5*rm*F(6,100))
    check("all_orders_series_denominator",5*rm*F(6,100)<1)
    check("complete_general_local_remainder",one+low+high<8*10**16,
          {"first":float(one),"orders_2_to_4":float(low),"all_higher":float(high)})
    check("scalar_local_remainder",F(103,100)*8*10**16<10**17)
    check("local_cauchy_factor",10**17*200**2==4*10**21)
    check("scalar_gluing_reserve",rm*6*10**11<7*10**11)
    check("gluing_cauchy_factor",7*10**11*200**2==28*10**15)

    J=F(10**8+1)
    qbound=4*J+3*(2*J)**2+4*2*dc*(2*J)+4*2*dc**2*(2*J)
    check("general_local_curvature_source",qbound<2*10**17,float(qbound))
    # |alpha|/C<100/(960*pi)<100/(960*3).
    actionbound=8*F(100,960*3)/(1-4*dc)*2*10**17
    check("common_contact_analytic_bound",actionbound<10**17,float(actionbound))
    eta=F(85,10**34)
    local=2*eta*F(2,3)*4*10**21
    gluing=2*eta*28*10**15*24*256**3
    contact=8*10**21*eta
    check("complete_general_interior_budget",local+gluing+contact<F(2,10**7),
          {"local":float(local),"actual_domain":float(gluing),
           "common_contact":float(contact),"total":float(local+gluing+contact)})


def numerical_quadratures():
    for name,fun,target in MOMENTUM_CASES:
        values=[]
        for nr,nt in ((64,48),(128,96)):
            nodes,weights=leggauss(nr)
            tt=(nodes+1)/2
            r=np.tan(np.pi*tt/2)
            jac=weights/2*(np.pi/2)*(1+r*r)
            angle=2*np.pi*np.arange(nt)/nt
            xx=r[:,None]*np.cos(angle)
            yy=r[:,None]*np.sin(angle)
            integrand=np.asarray(fun(xx,yy),dtype=complex)
            value=np.sum(jac*r*np.mean(integrand.real,axis=1))/(2*np.pi)
            values.append(float(value))
        tol=1e-10+1e-8*abs(target)
        check("independent_quadrature_"+name,
              all(abs(v-target)<tol for v in values) and abs(values[1]-values[0])<tol,
              {"target":target,"values":values})
        NUMERICS.append({"case":name,"target":target,"values":values})
    for i,fc in enumerate((.005+.004j,-.006+.003j,.001-.009j)):
        ratios=[]
        for M in (18.,23.):
            m=M*np.exp(fc)
            for r in (.0001,.002,.02,.1,1.):
                actual=max(abs(m*(kv(0,m*r)+kv(1,m*r))/(2*np.pi)),
                           abs(m*(kv(0,m*r)-kv(1,m*r))/(2*np.pi)))
                bound=1.02*np.sqrt(2)/np.pi*np.exp(-.98*M*r/np.sqrt(2))/r
                ratios.append(float(actual/bound))
        check(f"complex_mass_kernel_diagnostic_{i}",max(ratios)<1,
              {"phi_center":str(fc),"maximum_ratio":max(ratios)})


def finite_matrix_checks():
    mp.mp.dps=70
    M=mp.mpf(18)
    K=6*mp.matrix([[1,2,0,0],[-1,1,3,0],[0,-2,1,2],[0,0,-1,2]])
    D=mp.zeros(8)
    for i in range(4):
        for j in range(4):
            D[i,j+4]=K[i,j]
            D[j+4,i]=-K[i,j]
    base=[mp.mpf(v) for v in (".0002",".0003","-.0002",".0001")]*2
    u=mp.diag([mp.mpf(v) for v in (".2",".4","-.2",".3")]*2)
    w=mp.diag([mp.mpf(v) for v in ("-.1","-.3",".6",".2")]*2)
    u0,w0=u[0,0],w[0,0]
    ev=lambda A:A[0,0]+A[4,4]

    def data(aa=0,bb=0):
        vals=[base[i]+aa*u[i,i]+bb*w[i,i] for i in range(8)]
        mass=mp.diag([M*mp.exp(v) for v in vals])
        return (D+mass)**-1,mass,vals[0]

    def scalar(aa=0,bb=0):
        R,m,c=data(aa,bb)
        return -M*mp.exp(-c)*ev(R)

    def density(aa=0,bb=0):
        R,m,c=data(aa,bb)
        return -M*mp.exp(c)*ev(R)

    R,m,c=data()
    L=ev(R)
    Lu=-ev(R*m*u*R)
    Lw=-ev(R*m*w*R)
    Luw=ev(R*m*u*R*m*w*R+R*m*w*R*m*u*R-R*m*u*w*R)
    T=-M*mp.exp(-c)*L
    Tu=-M*mp.exp(-c)*(Lu-u0*L)
    Tw=-M*mp.exp(-c)*(Lw-w0*L)
    Tuw=-M*mp.exp(-c)*(Luw-u0*Lw-w0*Lu+u0*w0*L)
    Fuw=-M*mp.exp(c)*(Luw+u0*Lw+w0*Lu+u0*w0*L)
    restored=mp.exp(2*c)*(Tuw+2*u0*Tw+2*w0*Tu+4*u0*w0*T)
    check("scalar_density_chain_rule",abs(Fuw-restored)<mp.mpf("1e-60"))
    check("observation_values_are_nonzero",u0!=0 and w0!=0 and c!=0)
    for name,fun,target in (("scalar",scalar,Tuw),("density",density,Fuw)):
        errors=[]
        for hs in (".0001",".00005"):
            step=mp.mpf(hs)
            got=(fun(step,step)-fun(step,-step)-fun(-step,step)+fun(-step,-step))/(4*step**2)
            errors.append(abs(got-target))
        check(name+"_independent_mixed_difference",max(errors)<mp.mpf("1e-8"),
              {"target":mp.nstr(target,35),"errors":[mp.nstr(e,12) for e in errors]})
        check(name+"_difference_refinement",
              errors[1]<mp.mpf(".35")*errors[0] or max(errors)<mp.mpf("1e-45"))
        NUMERICS.append({"matrix_case":name,"target":mp.nstr(target,35),
                         "errors":[mp.nstr(e,12) for e in errors]})
    proper=mp.mpf(".00001")
    reference=-M*M*mp.besselk(0,M*proper)/mp.pi
    for fc in (mp.mpf("-.01"),mp.mpf(".007")):
        mass=M*mp.exp(fc)
        coordinate=proper*mp.exp(-fc)
        got=-M*mp.exp(-fc)*mass*mp.besselk(0,mass*coordinate)/mp.pi
        check("proper_distance_vacuum_"+str(fc),abs(got-reference)<mp.mpf("1e-60"))


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    initial=hashes()
    error=None
    try:
        exact_loop_checks()
        jet_and_action_checks()
        rational_bound_checks()
        numerical_quadratures()
        finite_matrix_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    final=hashes()
    ok=error is None and initial==final and all(row["passed"] for row in ROWS)
    record={
        "story":"S145","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
        "scipy":scipy.__version__,"input_sha256_before":initial,
        "input_sha256_after":final,"unchanged_inputs":initial==final,
        "checks":ROWS,"numerics":NUMERICS,"exception":error,
        "scope":"ordinary interior canonical massive scalar-source response "
                "for general conformal C5 directions; diagnostics are not continuum proof",
        "unresolved":["outer weighted response","global exterior-norm mapping",
                      "Weyl/LLL and shape/radius coupling","nonlinear exterior and physical gates"],
    }
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"]=(ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out=ROOT/"evidence"/f"{STEM}.json"
    history=json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs":[]}
    history["runs"].append(record)
    out.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    with (ROOT/"evidence"/f"{STEM}.log").open("a",encoding="utf-8",newline="\n") as log:
        log.write(json.dumps({"created_utc":record["created_utc"],"status":record["status"],
                              "inputs":initial})+"\n")
        for row in ROWS:
            log.write(json.dumps(row)+"\n")
        if error:
            log.write(error+"\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),
                      "failed":sum(not row["passed"] for row in ROWS)}))
    return 0 if ok else 1


if __name__=="__main__":
    sys.exit(main())
