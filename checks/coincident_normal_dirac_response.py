"""S144 exact local coefficients, quantitative bounds and independent diagnostics.

Continuum chart, resolvent and angular-limit arguments remain ordinary proofs.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import kv
import scipy

ROOT=Path(__file__).resolve().parents[1]
STEM="coincident-normal-dirac-response"
INPUTS=[
    "checks/coincident_normal_dirac_response.py",
    "notes/coincident-normal-dirac-response.md",
    "notes/curved-remote-dirac-response.md",
    "notes/flat-variable-dirac-three-point.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-curved-quantum-tt.md",
    "notes/stationary-reduced-parameter-form-bounds.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS=[]
NUMERICS=[]


def check(name, passed, detail=""):
    if any(row["name"]==name for row in ROWS):
        raise ValueError("duplicate check name")
    row={"name":name,"passed":bool(passed),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def exact_local_coefficients():
    x,y,a,b,c,B,alpha,sp,tp=s.symbols("x y a b c B alpha sp tp",real=True)
    gamma=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]])]
    N=s.eye(2)-s.I*(x*gamma[0]+y*gamma[1])
    den=1+x*x+y*y
    S=N/den
    H=s.Matrix([[a,c],[c,b]])
    hk=H*s.Matrix([x,y])
    khk=a*x*x+2*c*x*y+b*y*y
    P=-2*(a+b)*N*den+4*s.I*(hk[0]*gamma[0]+hk[1]*gamma[1])*den+8*khk*N
    differentiated=a*S.diff(x,2)+2*c*S.diff(x,y)+b*S.diff(y,2)
    check("momentum_second_derivative_identity",
          all(s.cancel(differentiated[i,j]-P[i,j]/den**3)==0
              for i in range(2) for j in range(2)))
    numerator=s.Poly(s.expand(-s.trace(P*N*P)/4),x,y)
    integral=s.Integer(0)
    for (ix,iy),coef in numerator.terms():
        if ix%2 or iy%2:
            continue
        i,j=ix//2,iy//2
        assert i+j<6
        moment=s.gamma(i+s.Rational(1,2))*s.gamma(j+s.Rational(1,2))*\
            s.gamma(6-i-j)/(4*s.pi**2*s.gamma(7))
        integral += coef*moment
    expected=(2*(a+b)**2+a*a+b*b+2*c*c)/(60*s.pi)
    check("independent_quadratic_born_tensor",
          s.simplify(integral-expected)==0,str(s.factor(integral)))
    check("traceless_hessian_cancellation",
          s.expand(-(B+(a+b)**2+2*(a*a+b*b+2*c*c))
                   +2*(2*(a+b)**2+a*a+b*b+2*c*c)
                   +B-3*(a+b)**2)==0)

    def bessel_square_moment(power,nu):
        ss=s.Integer(power+1)
        return 2**(ss-3)*s.gamma(ss/2)**2*s.gamma(ss/2+nu)*s.gamma(ss/2-nu)/s.gamma(ss)
    for p,target in ((3,-s.Rational(1,3)),(5,-s.Rational(8,15))):
        got=s.simplify(bessel_square_moment(p,0)-bessel_square_moment(p,1))
        check(f"bessel_difference_moment_{p}",got==target,str(got))
    lap=lambda f:s.diff(f,x,2)+s.diff(f,y,2)
    phi=(a*x*x+2*c*x*y+b*y*y)/2+B*x**4/24
    at0=lambda f:s.simplify(f.subs({x:0,y:0}))
    fourth=at0(lap(lap(s.exp(phi))))
    check("mass_exponential_fourth_jet",
          s.expand(fourth-(B+(a+b)**2+2*(a*a+b*b+2*c*c)))==0)
    u1,u2,w1,w2=s.symbols("u1 u2 w1 w2",real=True)
    total=phi+sp*(u1*x*x+u2*y*y)/2+tp*(w1*x*x+w2*y*y)/2
    # Euler-Lagrange derivative of 4 alpha exp(-2 phi)(Delta phi)^2.
    source=-8*alpha*s.exp(-2*total)*lap(total)**2+8*alpha*lap(s.exp(-2*total)*lap(total))
    source0=at0(source)
    check("direct_curvature_action_source",
          s.simplify(source0.subs({sp:0,tp:0})-8*alpha*(B-3*(a+b)**2))==0)
    mixed=s.diff(source0,sp,tp)
    check("restored_contact_mixed_variation",
          s.simplify(mixed+48*alpha*(u1+u2)*(w1+w2))==0)
    check("single_mode_action_source_sign",
          s.simplify((8*alpha*(B-3*(a+b)**2)).subs(alpha,-1/(960*s.pi))
                     +(B-3*(a+b)**2)/(120*s.pi))==0)


def exact_bounds():
    check("analytic_fifth_remainder",F(52,1000)*32**5 < 1750000)
    # Exponential Bell bounds with |zeta^(j)|<1 and exp(|zeta|)<1.001.
    bell=[0,1,2,5,15,52]
    reserve=[]
    for j in range(6):
        val=F(0)
        for i in range(j+1):
            seed=F(1026,1000) if i==0 else F(26,1000)*math.factorial(i)*64**i
            ex=F(1,1000) if j-i==0 else F(1001,1000)*bell[j-i]
            val += math.comb(j,i)*seed*ex
        reserve.append(val)
    check("finite_jet_product_reserve",max(reserve)<10**8,
          [float(v) for v in reserve])
    check("fourth_taylor_coefficient_reserve",
          F(26,1000)*32**4+reserve[4]/math.factorial(4)<10**6)
    check("inner_fifth_coefficient_cone",
          F(52,1000)*32**5+F(10**8,math.factorial(5))<10**12)
    check("outer_cubic_coefficient_cone",18*64<2*10**4)
    outer5=F(3,100)*64**5+9*64**3+2*10**4*64**2+10**6*64
    check("outer_fifth_coefficient_cone",outer5<10**12,str(outer5))
    check("cap_value_reserve",F(28,1000)/(1-F(28,1000))<F(3,100))
    # sqrt(2)<10/7, 1/sqrt(2)>7/10 give conservative rational norms.
    c_lower=F(7,10)-F(1,10)
    check("weighted_kernel_L1",2*F(10,7)/c_lower<5)
    check("weighted_kernel_L_three_halves_cubed",8*F(10,7)/(F(3,2)*c_lower)<13)
    check("two_edge_fifth_moment",F(8*24,1)/(3*F(7,5)**5)<16)
    one=16*10**12
    two=26*math.factorial(5)*10**5*(2*9*2*10**4)
    higher=F(26*math.factorial(6)*10**6*9**3*5,1)/(1-F(15,100))
    check("complete_local_remainder",F(one+two)+higher/18<2*10**14,
          {"one":one,"two":two,"higher_coefficient":float(higher)})
    check("local_two_parameter_cauchy",2*10**14*40**2==32*10**16)
    complex_gap=F(94,100)*(1-F(6,100)**2/2)-F(1,2)
    check("complex_weighted_resolvent",complex_gap>F(1,3))
    check("complex_flat_mass",1/(1-F(11,100))<F(113,100))
    check("complex_flat_mass_gradient",F(113,100)*(10+F(20,19)*F(6,100))<12)
    check("complex_squared_potential",F(113,100)**2+F(12,18)<3)
    # sqrt(pi)>sqrt(3)>17/10, density <= exp(.05)<20/19.
    check("improved_point_evaluation",F(30000,2)*F(10,17)*F(20,19)<10**4)
    check("improved_separated_row",3*10**4*F(18,17)<32*10**3)
    check("cutoff_gradient",F(15,8)*128*F(20,19)<256)
    check("cap_separation",F(95,100)/128>F(1,256))
    check("evaluation_disk_before_commutator",F(1,9*18)<F(95,100)/128)
    check("gluing_source_reserve",2*(32*10**3)**2*256<6*10**11)
    check("gluing_cauchy_reserve",6*10**11*40**2==96*10**13)
    E,d,C=s.symbols("E d C",positive=True)
    check("angular_inverse_cubic_tail",
          s.integrate(12*C/E**2,(E,18,s.oo))==2*C/3)
    check("angular_exponential_positive_majorant",
          s.integrate(4*C*d*E**3*s.exp(-d*E),(E,0,s.oo))==24*C/d**3)
    eta=F(85,10**34)
    remainder=2*eta*(F(2,3)*32*10**16+96*10**13*24*256**3)
    contact=14*eta
    check("normalized_remainder_budget",remainder<F(658,10**11),float(remainder))
    check("common_contact_budget",contact<F(12,10**32),float(contact))
    check("total_coincident_normal_budget",remainder+contact<F(7,10**9),
          {"remainder":float(remainder),"contact":float(contact)})


def numerical_born_tensor(H, nr, nt):
    # Independent direct matrix derivative, radial map k=tan(pi t/2).
    rr,rw=leggauss(nr)
    t=(rr+1)/2
    k=np.tan(np.pi*t/2)
    jac=(np.pi/2)*(1+k*k)*rw/2
    theta=2*np.pi*np.arange(nt)/nt
    kx=k[:,None]*np.cos(theta)
    ky=k[:,None]*np.sin(theta)
    I=np.eye(2,dtype=complex)
    gx=np.array([[0,1],[1,0]],dtype=complex)
    gy=np.array([[0,-1j],[1j,0]],dtype=complex)
    N=I-1j*(kx[...,None,None]*gx+ky[...,None,None]*gy)
    den=(1+kx*kx+ky*ky)[...,None,None]
    hkx=H[0,0]*kx+H[0,1]*ky
    hky=H[1,0]*kx+H[1,1]*ky
    khk=kx*hkx+ky*hky
    deriv=-2*np.trace(H)*N/den**2+\
        4j*(hkx[...,None,None]*gx+hky[...,None,None]*gy)/den**2+\
        8*khk[...,None,None]*N/den**3
    values=-np.trace(deriv@(N/den)@deriv,axis1=-2,axis2=-1).real/4
    return float(np.sum(jac*k*np.mean(values,axis=1))/(2*np.pi))


def numerical_checks():
    for p,target in ((3,-1/3),(5,-8/15)):
        val,err=quad(lambda r:r**p*(kv(0,r)**2-kv(1,r)**2),0,np.inf,
                     epsabs=1e-11,epsrel=1e-11)
        check(f"independent_bessel_quadrature_{p}",abs(val-target)<1e-9,
              {"value":val,"quadrature_estimate":err})
    for i,H in enumerate((np.diag([1.,0.]),np.array([[1.,.3],[.3,2.]]),
                          np.array([[1.,.7],[.7,-1.]]))):
        target=(2*np.trace(H)**2+np.trace(H@H))/(60*np.pi)
        values=[numerical_born_tensor(H,48,32),numerical_born_tensor(H,96,64)]
        tolerance=1e-10+1e-8*abs(target)
        check(f"independent_two_born_quadrature_{i}",
              all(abs(v-target)<tolerance for v in values)
              and abs(values[1]-values[0])<tolerance,
              {"H":H.tolist(),"target":target,"resolutions":values})
        NUMERICS.append({"H":H.tolist(),"target":target,"values":values})

    mp.mp.dps=70
    M=mp.mpf(18)
    K=6*mp.matrix([[1,2,0,0],[-1,1,3,0],[0,-2,1,2],[0,0,-1,2]])
    D=mp.zeros(8)
    for i in range(4):
        for j in range(4):
            D[i,j+4]=K[i,j]
            D[j+4,i]=-K[i,j]
    base=[mp.mpf(x) for x in ("0",".0003","-.0002",".0001")]*2
    u=mp.diag([mp.mpf(x) for x in ("0",".4","-.2",".3")]*2)
    w=mp.diag([mp.mpf(x) for x in ("0","-.3",".6",".2")]*2)
    ev=lambda A:A[0,0]+A[4,4]
    S=(D+M*mp.eye(8))**-1

    def operator(aa=0,bb=0):
        mass=mp.diag([M*mp.exp(base[i]+aa*u[i,i]+bb*w[i,i]) for i in range(8)])
        return (D+mass)**-1,mass

    def source(aa=0,bb=0):
        R,_=operator(aa,bb)
        return -M*ev(R-S)

    def first(bb=0):
        R,m=operator(bb=bb)
        return M*ev(R*m*u*R)

    def mixed(bb=0):
        R,m=operator(bb=bb)
        return M*ev(R*m*u*w*R-R*m*u*R*m*w*R-R*m*w*R*m*u*R)

    R,m=operator()
    h=m/M-mp.eye(8)
    term=S.copy()
    series=mp.zeros(8)
    for n in range(1,13):
        term=-M*term*h*S
        series+=term
    delta=max(abs(h[i,i]) for i in range(8))
    bound=2*delta**13/(1-delta)
    check("finite_matrix_complete_born_series",
          abs(-M*ev(series)-source())<bound,
          {"error":mp.nstr(abs(-M*ev(series)-source()),12),
           "geometric_tail_bound":mp.nstr(bound,12)})
    target=mixed()
    errors=[]
    for hs in (".0001",".00005"):
        step=mp.mpf(hs)
        got=(source(step,step)-source(step,-step)-source(-step,step)+source(-step,-step))/(4*step**2)
        errors.append(abs(got-target))
    check("independent_coincident_source_difference",max(errors)<mp.mpf("1e-8"),
          {"target":mp.nstr(target,35),"errors":[mp.nstr(e,12) for e in errors]})
    check("coincident_difference_refinement",
          errors[1]<mp.mpf(".35")*errors[0] or max(errors)<mp.mpf("1e-45"))
    endpoint=mp.mpf(".001")
    integrated=mp.quad(mixed,[0,endpoint])
    change=first(endpoint)-first()
    check("finite_background_source_derivative_change",
          abs(change-integrated)<mp.mpf("1e-30"),
          {"change":mp.nstr(change,35),"quadrature":mp.nstr(integrated,35),
           "error":mp.nstr(abs(change-integrated),12)})
    NUMERICS.append({"matrix_mixed_source":mp.nstr(target,35),
                     "central_errors":[mp.nstr(e,12) for e in errors],
                     "finite_path_error":mp.nstr(abs(change-integrated),12)})


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    initial=hashes()
    error=None
    try:
        exact_local_coefficients()
        exact_bounds()
        numerical_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    final=hashes()
    ok=error is None and initial==final and all(r["passed"] for r in ROWS)
    record={
        "story":"S144","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
        "scipy":scipy.__version__,"input_sha256_before":initial,
        "input_sha256_after":final,"unchanged_inputs":initial==final,
        "checks":ROWS,"numerics":NUMERICS,"exception":error,
        "scope":"ordinary coincident massive point-source proof on the fixed "
                "interior normal-coordinate sector; numerical tests are diagnostics",
        "unresolved":["general value/first-jet coordinate return",
                      "outer finite-background and exterior-norm transfer",
                      "Weyl/LLL and coupled exterior","physical completion gates"],
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
                      "failed":sum(not r["passed"] for r in ROWS)}))
    return 0 if ok else 1


if __name__=="__main__":
    sys.exit(main())
