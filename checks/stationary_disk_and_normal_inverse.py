"""Finite identities for the global disk reduction and coupled outer normal inverse.

Heat localization, Mellin continuation, digamma estimates and weighted
multiplier mapping remain ordinary proofs. No global inverse is tested.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/stationary-disk-and-normal-inverse.json"
LOG = ROOT / "evidence/stationary-disk-and-normal-inverse.log"
INPUTS = [
    "checks/stationary_disk_and_normal_inverse.py",
    "notes/stationary-disk-and-normal-inverse.md",
    "notes/stationary-shape-logarithmic-response.md",
    "notes/stationary-weighted-conformal-range.md",
    "notes/stationary-fast-phase-variation.md",
    "notes/stationary-exterior-boundary-weyl.md",
    "notes/spherical-anomaly-action.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-same-source-boundary-charge.md",
    "notes/uniform-dirac-reference.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-paired-dirac-phase.md",
]
CONTRACT = {
    "acceptance": "every symbolic residual exactly zero and rational inequality true",
    "numerical_tolerances": "none",
    "scope": "finite identities for global disk equations and coupled reflecting outer normal inverse",
    "ordinary_arguments": ["global smooth disk and bolt parametrix", "even-double heat identity and normal localization",
        "Mellin continuation and digamma-series bounds", "weighted normal multiplier inverse"],
    "unresolved": ["actual-minus-normal maps", "mass extension", "global kernel/range"],
    "excluded": ["exact RN solution", "physical instability or nonexistence",
        "full coupled index or inverse", "actual compatibility Jacobian",
        "Planck frequency numerical shooting", "physical gate closure"],
    "failure_policy": "append every run and retain failed input bytes",
    "hash_policy": "SHA256 of raw UTF8/LF/noBOM bytes, no normalization",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row={"label":label,"passed":bool(condition),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row,sort_keys=True),flush=True)


def zero(label,expression):
    value=S.factor(S.together(S.expand(expression)))
    check(label,value==0,value)


def disk_and_ward():
    x,B,kappa,lam=S.symbols("x B kappa lam",positive=True)
    phi=S.Function("phi")(x);v=S.Function("v")(x)
    radius=B*S.exp(phi)
    F=B**2*S.exp(-2*v)*S.diff(phi,x)**2
    alpha=S.log(kappa*x/B)+phi+2*v-S.log(S.diff(phi,x))
    zero("disk_areal_radial_metric",S.diff(radius,x)**2/F-S.exp(2*phi+2*v))
    zero("disk_time_metric",S.exp(2*alpha)*F/kappa**2-S.exp(2*phi+2*v)*x**2)
    zero("disk_normalization_ODE",kappa*S.exp(-alpha)*S.diff(radius,x)/F-1/x)
    ph,ph2,vh=S.symbols("ph ph2 vh",real=True)
    rh=B*S.exp(ph)
    # Leading coefficients in xi, before the singular areal quotient.
    Flead=4*B**2*S.exp(-2*vh)*ph2**2
    rlead=rh*ph2
    ah=kappa*S.exp(ph+2*vh)/(2*B*ph2)
    zero("bolt_temperature",ah*(Flead/rlead)/2-kappa)
    a=S.Function("a")(x);r=S.Function("r")(x);e=S.Function("e")(x)
    ward=S.diff(e,x)+2*S.diff(a,x)*e/a+2*S.diff(r,x)*e/r
    zero("orbit_trace_Ward_integrating_factor",S.diff(a**2*r**2*e,x)-a**2*r**2*ward)
    const=S.symbols("const")
    zero("orbit_trace_residual_solution",ward.subs(e,const/(a**2*r**2)).doit())
    dh,df,dfp,fpp=S.symbols("dh df dfp fpp")
    zero("moving_bolt_pullback_value",(df+2*kappa*dh).subs(df,-2*kappa*dh))
    zero("moving_bolt_pullback_derivative",(dfp+fpp*dh).subs(dfp,-fpp*dh))
    m=S.symbols("m")
    z=lam*x/kappa+lam**2*x**3/(3*kappa**3)-m*lam**3*x**4/(2*kappa**4)
    original=kappa*(z/lam-z**3/(3*lam**2)+m*z**4/(2*lam**2))
    zero("outer_coordinate_series_reversion",S.series(original-x,x,0,5).removeO())
    phi_correction=S.series(-S.log(z/(lam*x/kappa)),x,0,4).removeO()
    zero("outer_radius_log_series",phi_correction+lam*x**2/(3*kappa**2)-m*lam**2*x**3/(2*kappa**3))
    vc=S.series(x+S.log(1+z**2/lam-2*m*z**3/lam)/2,x,0,4).removeO()
    zero("outer_orbit_series",vc-x-lam*x**2/(2*kappa**2)+m*lam**2*x**3/kappa**3)
    zero("one_correlated_mass_cubic",-lam**2/kappa**3+2*lam**2/(2*kappa**3))
    sig,tau=S.symbols("sig tau")
    T2=(2*sig+2*tau)**2+2*(2*sig)**2
    zero("disk_transverse_symbol",T2-12*sig**2-8*sig*tau-4*tau**2)
    ell=S.symbols("ell")
    poly=T2+8*tau**2*(ell-R(23,15))
    zero("disk_mixed_square",poly-12*(sig+tau/3)**2-8*tau**2*(ell-R(6,5)))


def clifford_and_curvature():
    pa=[S.Matrix([[0,1],[1,0]]),S.Matrix([[0,-S.I],[S.I,0]]),S.diag(1,-1)]
    ident=S.eye(4)
    gamma=[S.kronecker_product(pa[2],S.eye(2))]+[S.kronecker_product(pa[0],p) for p in pa]
    chi=gamma[0]
    chir=gamma[0]*gamma[1]*gamma[2]*gamma[3]
    check("four_dimensional_Clifford_algebra",all(gamma[i]*gamma[j]+gamma[j]*gamma[i]==(2*ident if i==j else S.zeros(4)) for i in range(4) for j in range(4)))
    check("normal_projector_reflection",all(chi*(gamma[0]*gamma[i])*chi==-gamma[0]*gamma[i] for i in range(1,4)))
    check("chirality_commutes_even_square",all(chir*(gamma[0]*gamma[i])==(gamma[0]*gamma[i])*chir for i in range(1,4)))
    check("chirality_reverses_image_trace",chir*chi==-chi*chir and chir**2==ident)
    p1,p2,q1,q2,t,w1,w2,h=S.symbols("p1 p2 q1 q2 t w1 w2 h")
    pp=[p1,p2,-p1-p2];qq=[q1,q2,-q1-q2]
    kn=[-(t*w2+h*qq[i]+t*w1*h*pp[i]+h**2*pp[i]**2) for i in range(3)]
    kt=lambda i,j:-(t*w1+h*pp[i])*(t*w1+h*pp[j])
    euler=8*(kn[0]*kt(1,2)+kn[1]*kt(0,2)+kn[2]*kt(0,1))
    U=sum(p*p for p in pp);Up=2*sum(p*q for p,q in zip(pp,qq))
    zero("Euler_TT_quadratic_density",S.expand(euler).coeff(h,2)+4*t*(w2*U+w1*Up))
    zero("Euler_TT_linear_density",S.expand(euler).coeff(h,1))
    scalar=-h**2*U-6*t*w2-6*t**2*w1**2
    R2=S.expand(scalar**2-(h**2*U)**2)
    zero("R_squared_TT_quadratic_density",R2.coeff(h,2)-12*(t*w2+t**2*w1**2)*U)
    bp=S.symbols("bp")
    zero("TT_beta_normalization",(bp/3)**2+2*(-bp/6)**2-bp**2/6)
    G,lam=S.symbols("G lam",positive=True)
    aE=R(11,5760)/S.pi**2;d=R(1,480)/S.pi**2
    AT=1/(48*S.pi*G*lam)+2*(aE+d)/3
    zero("Euler_R2_Einstein_TT_coefficient",AT-1/(48*S.pi*G*lam)-R(23,8640)/S.pi**2)
    zero("dimensionless_TT_Einstein_ratio",AT/d-10*S.pi/(G*lam)-R(23,18))


def normal_algebra():
    s,G,lam,H,omega=S.symbols("s G lam H omega")
    z=S.symbols("z",positive=True)
    q=s*(s-3);P4=s*(s-1)*(s-2)*(s-3)
    zero("normal_quartic_factor",P4-q*(q+2))
    zero("Mellin_adjoint_polynomial",P4-P4.subs(s,3-s))
    eps=G/(60*S.pi);aa=1+R(11,3)*eps*lam
    d0=1/(480*S.pi**2)
    AT=1/(48*S.pi*G*lam)+R(23,8640)/S.pi**2
    NT=d0*P4*H-AT*q
    JC=3*(aa+eps*lam*q)*(q-4)/(4*S.pi*G*lam)
    zero("scalar_oscillatory_root",(aa+eps*lam*q).subs(aa,eps*lam*(omega**2+R(9,4))).subs(s,R(3,2)+S.I*omega))
    zero("scalar_growing_root",JC.subs(s,-1))
    zero("scalar_fast_root",JC.subs(s,4))
    eta=z**(s+1)/(3*(s+1))
    zero("normal_radial_gauge_primitive",S.powsimp(S.diff(eta,z))-z**s/3)
    zero("disk_to_TT_conformal_coefficient",R(1,3)+2/(3*(s+1))-(s+3)/(3*(s+1)))
    C=S.Matrix([[1,(s+3)/(3*(s+1))],[0,2]])
    adj=C.subs(s,3-s).T
    disk=adj*S.diag(JC,NT)*C
    zero("reflected_normal_determinant",S.det(disk)-4*JC*NT)
    zero("mass_cubic_is_pure_TT",1+(3+3)/(S.Integer(3)*(3+1))*(-2))
    old=(s-3)/(6*(s+1))
    old_adj=old.subs(s,3-s)
    zero("Einstein_shape_diagonal_cancellation",NT+old_adj*JC*old-d0*(P4*H+q*(q-4)/6))
    zero("TT_to_original_shape_constant",-R(6,5)+R(1,6)+R(31,30))
    zero("scalar_contour_distance_gap",R(1,2)*R(9,2)*2*2/(80*S.pi**2)-R(9,80)/S.pi**2)
    zero("mass_root_normal",NT.subs(s,3))
    zero("source_root_normal",NT.subs(s,0))
    zero("flat_mass_log_residue",S.limit(-P4/(s-3),s,3)+6)
    zero("Weyl_mass_contact",S.diff(S.log(z)*S.diff(z**3,z,2),z,2)-6/z)
    zero("first_pole_zero_product",S.limit(-P4/(s-1),s,1)+2)
    zero("second_pole_zero_product",S.limit(P4/(s-2),s,2)+2)
    u=S.symbols("u")
    fourier=-2*S.gamma(u+1)*S.sin(S.pi*u/2)
    logder=S.diff(fourier,u)/fourier-S.polygamma(0,u+1)-S.pi*S.cos(S.pi*u/2)/(2*S.sin(S.pi*u/2))
    zero("Fourier_power_log_derivative",logder)
    zero("Weyl_log_cancellation",(H-S.log(z))-(S.log(S.Symbol("length"))-S.log(z))-(H-S.log(S.Symbol("length"))))


def contour_guards():
    t=S.symbols("t",real=True)
    st=1/(t**2+R(1,4))+1/(t**2+R(9,4))
    q=(R(7,2)+S.I*t)*(R(1,2)+S.I*t)
    zero("contour_quadratic",S.expand(q)-(R(7,4)-t**2+4*S.I*t))
    h=S.symbols("h",real=True)
    product=(q+2)*(h-S.I*t*st)
    zero("contour_imaginary_identity",S.im(S.expand(product))-t*(4*h-(R(15,4)-t**2)*st))
    zero("contour_real_identity",S.re(S.expand(product))-(R(15,4)-t**2)*h-4*t**2*st)
    check("benchmark_scale_lower_integer_guard",3**18<R(2000**3,4))
    zero("scheme_constant_lower",9-R(6,5)-R(39,5))
    zero("scheme_constant_upper",12-R(6,5)-R(54,5))
    zero("digamma_denominator_bound",4+R(4,9)-R(40,9))
    zero("positive_imaginary_gap",4*R(29,5)-R(15,4)*R(40,9)-R(98,15))
    check("imaginary_gap_positive",R(98,15)>0)
    zero("small_frequency_real_guard",R(15,4)*R(84,5)+8-71)
    check("real_guard_below72",71<72)
    check("Einstein_ratio_gap",30*2000**10>144)
    check("real_gap_positive",144-72>0)
    check("all_scalar_root_distances_positive",all(v>0 for v in [R(1,2),R(9,2),2,2]))


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before={name:(ROOT/name).read_bytes() for name in INPUTS}
    existing=json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run={"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
         "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
         "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run,sort_keys=True),flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("inputs_strict_UTF8_LF",all(raw.decode("utf-8") is not None and b"\r" not in raw
                  and not raw.startswith(b"\xef\xbb\xbf") for raw in before.values()))
            disk_and_ward()
            clifford_and_curvature()
            normal_algebra()
            contour_guards()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"]=hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"],sort_keys=True),flush=True)
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"]=ROWS
        run["failed_checks"]=[row for row in ROWS if not row["passed"]]
        code=int(bool(run["failed_checks"]))
        run["status"]="failed" if code else "passed"
        run["terminal_exit_status"]=code
        run["finished_utc"]=datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"]={name:base64.b64encode(raw).decode("ascii") for name,raw in before.items()}
            for name,raw in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(raw).hexdigest(),
                                  "base64":base64.b64encode(raw).decode("ascii")},sort_keys=True))
        print("TERMINAL "+json.dumps({"status":run["status"],"checks":len(ROWS),
              "failures":len(run["failed_checks"]),"terminal_exit_status":code,
              "finished_utc":run["finished_utc"]},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),
                      "receipt":str(RECEIPT)},sort_keys=True))
    return code


if __name__=="__main__":
    sys.exit(main())
