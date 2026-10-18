"""Finite identities for exact compact shape response and mixed interior ellipticity.

Heat localization, pseudodifferential estimates and endpoint/global
solvability remain ordinary arguments. No physical stability is tested.
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
RECEIPT = ROOT / "evidence/stationary-shape-logarithmic-response.json"
LOG = ROOT / "evidence/stationary-shape-logarithmic-response.log"
INPUTS = [
    "checks/stationary_shape_logarithmic_response.py",
    "notes/stationary-shape-logarithmic-response.md",
    "notes/stationary-weighted-conformal-range.md",
    "notes/stationary-fast-phase-variation.md",
    "notes/stationary-coupled-asymptotic-obstruction.md",
    "notes/spherical-anomaly-action.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-paired-dirac-phase.md",
]
CONTRACT = {
    "acceptance": "every symbolic residual exactly zero and rational inequality true",
    "numerical_tolerances": "none",
    "scope": "finite identities for compact shape heat Hessian and coupled interior mixed-log symbol",
    "ordinary_arguments": ["compact heat Hessian and common Mellin continuation", "interior localization with retained global contacts",
        "mixed-log parametrix and Garding estimate", "packet function-space obstruction",
        "global shape/constraint maps and moving endpoint domains"],
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


def compact_vertices():
    u,B,f=S.symbols("u B f",positive=True)
    fp,v1,v2,d1,d2,t1,t2=S.symbols("fp v1 v2 d1 d2 t1 t2",real=True)
    ft=f+t1*v1+t2*v2
    fpt=fp+t1*d1+t2*d2
    A=u*S.sqrt(ft)/B
    H=u/(B*S.sqrt(ft))
    c=fpt/(4*ft)-1/(2*u)
    at=lambda z:z.subs({t1:0,t2:0})
    first=lambda z:at(S.diff(z,t1))
    mixed=lambda z:at(S.diff(z,t1,t2))
    a0,h0,c0=at(A),at(H),at(c)
    b1,b2=v1/f,v2/f
    bp1,bp2=d1/f-v1*fp/f**2,d2/f-v2*fp/f**2
    prodprime=bp1*b2+b1*bp2
    zero("product_orbit_density",(B**2*ft/u**2)*(B**2/(u**2*ft))-B**4/u**4)
    zero("fixed_volume_first_variation",S.diff((B**2*ft/u**2)*(B**2/(u**2*ft)),t1))
    zero("Dirac_radial_first",first(A)-a0*b1/2)
    zero("Dirac_time_first",first(H)+h0*b1/2)
    zero("spin_connection_first",first(c)-bp1/4)
    zero("Dirac_radial_contact",mixed(A)+a0*b1*b2/4)
    zero("Dirac_time_contact",mixed(H)-3*h0*b1*b2/4)
    zero("spin_connection_contact",mixed(c)+prodprime/4)
    zero("full_radial_zero_order_first",first(A*c)-a0*(b1*c0/2+bp1/4))
    zero("full_radial_zero_order_contact",mixed(A*c)-a0*(-b1*b2*c0/4-prodprime/8))
    D,V1,V2,V12,Ad,J1,J2,J12=S.symbols("D V1 V2 V12 Ad J1 J2 J12",commutative=False)
    square=S.expand((Ad+t1*J1+t2*J2+t1*t2*J12)*(D+t1*V1+t2*V2+t1*t2*V12))
    zero("positive_square_first_vertex",S.expand(square).coeff(t1,1).coeff(t2,0)-(J1*D+Ad*V1))
    zero("positive_square_mixed_contact",S.expand(square).coeff(t1,1).coeff(t2,1)-(J12*D+Ad*V12+J1*V2+J2*V1))
    q,x,y,z,t=S.symbols("q x y z t",positive=True)
    heat=S.exp(-t*(q+t1*x+t2*y+t1*t2*z))
    zero("scalar_Duhamel_time_ordering",at(S.diff(heat,t1,t2))-S.exp(-t*q)*(t**2*x*y-t*z))


def principal_hessian():
    sig,b,p,ell,z=S.symbols("sig b p ell z",real=True)
    p=S.symbols("p",positive=True)
    pi=S.pi
    K=S.diag(2*sig+b,2*sig-b,2*sig,2*sig)
    proj=S.diag(1,0,1,1)
    T1=S.trace(K*proj)
    T2=S.trace(K*proj*K*proj)
    phi=proj*K*proj-T1*proj/3
    zero("transverse_trace",T1-6*sig-b)
    zero("transverse_quadratic_trace",T2-12*sig**2-4*sig*b-b**2)
    zero("tracefree_projection_logarithm",3*T2-T1**2-2*b**2)
    zero("shape_projection_norm",S.trace(K*phi)-2*b**2/3)
    zero("pure_conformal_logarithm_vanishes",S.trace(K*phi).subs(b,0))
    zero("principal_scalar_curvature",p**2*(S.trace(K)-K[1,1])-p**2*(6*sig+b))
    eta=2*sum(R(1,j) for j in range(1,6))-sum(R(1,j) for j in range(1,3))
    zero("even_dimension_eta",eta-R(46,15))
    W=p**(4-2*z)/(16*pi**2)*((2*z-3)*T2+T1**2)
    quadratic=(S.diff(W,z).subs(z,0)+eta*W.subs(z,0))/120
    expected=p**4/(960*pi**2)*(12*sig**2+4*sig*b+b**2+2*b**2*(S.log(p)-eta/2))
    zero("zeta_to_fermionic_Hessian",quadratic-expected)
    leading=S.expand(expected).coeff(S.log(p))
    zero("universal_shape_log_coefficient",leading-p**4*b**2/(480*pi**2))
    G=S.symbols("G",positive=True)
    eps=G/(60*pi)
    hss=S.diff(expected,sig,2)/2
    hsb=S.diff(expected,sig,b)/2
    hbb=S.diff(expected,b,2)/2
    zero("Euclidean_conformal_trace_sign",hss-6*eps*p**4/(8*pi*G))
    zero("Euclidean_shape_trace_cross",hsb-eps*p**4/(8*pi*G))
    zero("conformal_principal_coefficient",hss-p**4/(80*pi**2))
    zero("cross_matrix_coefficient",hsb-p**4/(480*pi**2))
    theta=S.symbols("theta",real=True)
    target=p**4*theta**2/(80*pi**2)+p**4*b**2/(480*pi**2)*(S.log(p)+R(1,3)-eta/2)
    zero("mixed_square_completion",expected.subs(sig,theta-b/6)-target)
    zero("ordinary_shape_constant",R(1,3)-eta/2+R(6,5))
    zero("conformal_Schur_retains_shape_log",hbb-hsb**2/hss-p**4/(480*pi**2)*(S.log(p)-R(6,5)))
    zero("curvature_square_scalar_direction",(6*(theta-b/6)+b)**2-36*theta**2)
    J=S.Matrix([[1,-R(1,6)],[0,1]])
    mat=S.Matrix([[hss,hsb],[hsb,hbb]])
    diag=J.T*mat*J
    zero("transformed_matrix_cross",diag[0,1])
    zero("transformed_matrix_symmetry",diag[1,0])
    fs,fb=S.symbols("fs fb")
    outputs=J.T*S.Matrix([fs,fb])
    zero("transpose_output_conformal",outputs[0]-fs)
    zero("transpose_output_shape",outputs[1]-fb+fs/6)
    mu=S.symbols("mu",positive=True)
    scaled=expected.subs(S.log(p),S.log(p)-S.log(mu))
    zero("scale_shift_only_shape_order_four",scaled-expected+p**4*b**2*S.log(mu)/(480*pi**2))
    check("positive_conformal_symbol",S.ask(S.Q.positive(R(1,80)/pi**2)) is True)
    check("positive_shape_log_symbol",S.ask(S.Q.positive(R(1,480)/pi**2)) is True)


def mixed_inverse_and_energy():
    p,d,a,L,c,delta=S.symbols("p d a L c delta",positive=True)
    theta,b=S.symbols("theta b",real=True)
    leading=S.diag(a*p**4,d*p**4*(L+c))
    inverse=S.diag(1/(a*p**4),1/(d*p**4*(L+c)))
    product=leading*inverse
    zero("reciprocal_conformal_symbol",product[0,0]-1)
    zero("reciprocal_log_shape_symbol",product[1,1]-1)
    square=(S.sqrt(delta)*p**2*theta-p*L*b/(2*S.sqrt(delta)))**2
    zero("mixed_remainder_Young_identity",delta*p**4*theta**2+p**2*L**2*b**2/(4*delta)-p**3*L*theta*b-square)
    zero("shape_remainder_energy_ratio",(p**2*L**2)/(p**4*L)-L/p**2)
    zero("ordinary_shape_energy_ratio",p**4/(p**4*L)-1/L)


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
            compact_vertices()
            principal_hessian()
            mixed_inverse_and_energy()
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
