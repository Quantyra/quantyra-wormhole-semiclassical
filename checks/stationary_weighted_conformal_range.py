"""Exact identities for the weighted conformal trace range and coupled reduction.

Endpoint/Fredholm and domain arguments remain ordinary proofs. No
nonlocal shape inverse, actual Jacobian or physical matching is tested.
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
RECEIPT = ROOT / "evidence/stationary-weighted-conformal-range.json"
LOG = ROOT / "evidence/stationary-weighted-conformal-range.log"
INPUTS = [
    "checks/stationary_weighted_conformal_range.py",
    "notes/stationary-weighted-conformal-range.md",
    "notes/stationary-trace-range-condition.md",
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
    "scope": "finite identities for the off-shell weighted conformal trace range and constrained coupled reduction",
    "ordinary_arguments": ["bolt Volterra and Holder regularity", "weighted end estimates and matching",
        "Fredholm index and adjoint range sufficiency", "physical logarithmic bootstrap",
        "full state/domain and nonlocal shape mapping"],
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


def conformal_coefficients():
    eps,lam,lam0,a=S.symbols("eps lam lam0 a",positive=True)
    scalar,boxR,Weyl,J2,mag=S.symbols("scalar boxR Weyl J2 mag")
    sig,bs,bbs,gr,jh=S.symbols("sig bs bbs gr jh")
    dR=-6*bs-2*scalar*sig
    dBoxR=-6*bbs-2*scalar*bs-2*gr-4*boxR*sig
    E=scalar+12*lam0+eps*(boxR-R(11,72)*scalar**2+R(7,12)*Weyl+R(11,6)*J2-20*mag)
    Z=scalar+12*lam
    centered=(1+R(11,3)*eps*lam)*Z+eps*boxR-30*eps*(R(2,3)*mag-R(7,360)*Weyl-R(11,180)*J2+R(11,2160)*Z**2)
    zero("trace_vacuum_centering",E.subs(lam0,lam+R(11,6)*eps*lam**2)-centered)
    dE=dR+eps*(dBoxR-R(11,36)*scalar*dR+R(7,12)*(-4*sig*Weyl)
               +R(11,6)*(-4*sig*J2-4*jh)-20*(-4*sig*mag))
    V=-2*scalar+R(11,18)*eps*scalar**2-4*eps*boxR+80*eps*mag-R(7,3)*eps*Weyl-R(22,3)*eps*J2
    target=-6*eps*bbs+(-6-eps*scalar/6)*bs-R(22,3)*eps*jh-2*eps*gr+V*sig
    zero("exact_conformal_divergence_operator",dE-target)
    zero("potential_retains_off_shell_defect",V-(2*scalar+48*lam0-4*E))
    zero("divergence_tensor_coefficient",-eps/6-R(22,3)*eps/4+2*eps)
    j1,j2,j3,h1,h2,h3,h4=S.symbols("j1 j2 j3 h1 h2 h3 h4")
    js=[j1,j2,j3,-j1-j2-j3]
    hs=[h1,h2,h3,h4]
    contraction=sum(j*(-2*h+sum(hs)/2) for j,h in zip(js,hs))
    zero("traceless_Ricci_square_variation",2*contraction+4*sum(j*h for j,h in zip(js,hs)))
    t=S.symbols("t")
    zero("conformal_magnetic_invariant",S.diff(S.exp(-4*t*sig)*mag,t).subs(t,0)+4*sig*mag)
    dtrace=S.symbols("dtrace")
    zero("trace_density_variation",S.diff(S.exp(4*t*sig)*(E+t*dtrace),t).subs(t,0)-dtrace-4*sig*E)
    zero("AdS_second_order_coefficient",(-6-eps*scalar/6).subs(scalar,-12*lam)-(24*eps*lam-6*(1+R(11,3)*eps*lam)))
    zero("AdS_potential",(2*scalar+48*lam0).subs({scalar:-12*lam,lam0:lam+R(11,6)*eps*lam**2})-24*lam*(1+R(11,3)*eps*lam))


def radial_and_indicial():
    r,eps,lam,a=S.symbols("r eps lam a",positive=True)
    w,p,b,V=[S.Function(n)(r) for n in ["w","p","b","V"]]
    u,v=[S.Function(n)(r) for n in ["u","v"]]
    box=lambda z:S.diff(p*S.diff(z,r),r)/w
    op=lambda z:-6*eps*box(box(z))+S.diff(p*b*S.diff(z,r),r)/w+V*z
    concom=-6*eps*p*(u*S.diff(box(v),r)-S.diff(u,r)*box(v)+box(u)*S.diff(v,r)-S.diff(box(u),r)*v)+p*b*(u*S.diff(v,r)-S.diff(u,r)*v)
    zero("radial_Green_identity",S.diff(concom,r)-w*(u*op(v)-v*op(u)))
    zero("radial_second_order_product",S.diff(p*b*S.diff(v,r),r)/w-b*box(v)-(p/w)*S.diff(b,r)*S.diff(v,r))
    forcing=S.Function("forcing")(r)
    volterra_rhs=(b*box(v)+(p/w)*S.diff(b,r)*S.diff(v,r)+V*v-forcing)/(6*eps)
    zero("bolt_Volterra_equation_sign",-6*eps*(box(box(v))-volterra_rhs)-(op(v)-forcing))
    n,c,y=S.symbols("n c y",positive=True)
    near=lambda z:c*S.diff(y*S.diff(z,y),y)
    zero("bolt_principal_indicial_polynomial",near(near(y**n))-c**2*n**2*(n-1)**2*y**(n-2))
    s,om=S.symbols("s om",real=True)
    ind=-6*lam*(a+eps*lam*s*(s+3))*(s*(s+3)-4)
    zero("growing_indicial_root",ind.subs(s,1))
    zero("fast_indicial_root",ind.subs(s,-4))
    zero("oscillatory_indicial_root",ind.subs({a:eps*lam*(om**2+R(9,4)),s:-R(3,2)+S.I*om}))
    check("noncritical_weight_three",-4+3<0 and 1+3>0 and -R(3,2)+3>0)
    zero("resonant_derivative",S.diff(ind,s).subs(s,-4)-30*lam*(a+4*eps*lam))
    check("resonance_coefficient_positive",S.ask(S.Q.positive(30*lam*(a+4*eps*lam))) is True)
    limbox=lambda z:lam*(r**2*S.diff(z,r,2)+4*r*S.diff(z,r))
    limop=lambda z:-6*(a*(limbox(z)-4*lam*z)+eps*limbox(limbox(z)-4*lam*z))
    zero("indicial_operator",limop(r**s)-ind*r**s)
    zero("one_log_resonance",limop(r**(-4)*S.log(r))-30*lam*(a+4*eps*lam)*r**(-4))
    def green(x,z):
        return -6*eps*lam*r**4*(x*S.diff(limbox(z),r)-S.diff(x,r)*limbox(z)+limbox(x)*S.diff(z,r)-S.diff(limbox(x),r)*z)+lam*r**4*(24*eps*lam-6*a)*(x*S.diff(z,r)-S.diff(x,r)*z)
    zero("fast_growing_nonzero_pairing",green(r,r**(-4))-30*lam*(a+4*eps*lam))
    zero("oscillatory_fast_pairing",green(r**s,r**(-4)).subs(a,-eps*lam*s*(s+3)))
    zero("fast_self_pairing",green(r**(-4),r**(-4)))
    check("moment_integrable_weight_three",2-R(3,2)-3<-1)


def coordinate_constraint_blocks():
    u,r,eps,lam,a=S.symbols("u r eps lam a",positive=True)
    sig=S.Function("sig")(u)
    f=S.Function("f")(u)
    radial=u*S.exp(sig)
    q=1+u*S.diff(sig,u)
    F=f*q**2
    alpha=sig-S.log(q)
    zero("conformal_radius_derivative",S.diff(radial,u)-S.exp(sig)*q)
    zero("conformal_time_metric",S.exp(2*alpha)*F-S.exp(2*sig)*f)
    zero("conformal_radial_metric",S.diff(radial,u)**2/F-S.exp(2*sig)/f)
    zero("integral_coordinate_ODE",1/S.diff(radial,u)-S.exp(alpha)*u**2/radial**2)
    ff,fp,al,ap,uu,up=S.symbols("ff fp al ap uu up")
    fdot=S.exp(2*al)*uu**2/r**2*(fp+2*ff*(ap+up/uu-1/r))
    fu=S.cancel(fdot/up).subs(up,S.exp(al)*uu**2/r**2)
    zero("temperature_preservation",fu.subs(ff,0)-S.exp(al)*fp)
    kappa,duh,dfh,dfph,fpp=S.symbols("kappa duh dfh dfph fpp")
    zero("moving_bolt_value",(dfh+2*kappa*duh).subs(dfh,-2*kappa*duh))
    zero("moving_bolt_temperature",(dfph+fpp*duh).subs(dfph,-fpp*duh))
    e=S.Function("e")(r)
    bp=S.symbols("bp")
    ward=S.diff(e,r)+bp*(e-e)+2/r*(e-(-e))
    zero("Ward_trace_anisotropy_residual",ward-S.diff(e,r)-4*e/r)
    const=S.symbols("const")
    zero("surviving_radial_constraint_mode",ward.subs(e,const/r**4).doit())
    zero("radial_constraint_identifies_constant",r**4*(const/r**4)-const)
    pu=S.Function("pu")(r)
    F0=1+lam*r**2
    boxpu=F0*S.diff(pu,r,2)+(S.diff(F0,r)+2*F0/r)*S.diff(pu,r)
    Dr=F0*S.diff(pu,r,2)+S.diff(F0,r)*S.diff(pu,r)/2-boxpu+3*lam*pu
    Dt=S.diff(F0,r)*S.diff(pu,r)/2-boxpu+3*lam*pu
    zero("AdS_shape_cross_operator",-2*(Dr-Dt)/F0+2*S.diff(pu,r,2))
    s=S.symbols("s")
    ps=a+eps*lam*s*(s+3)
    cross=-2*s*(s-1)*ps
    zero("fast_trace_root_shape_residual",cross.subs(s,-4)+40*(a+4*eps*lam))
    check("shape_cross_symbol_nonzero",S.ask(S.Q.positive(40*(a+4*eps*lam))) is True)
    C,D,G,T,B,v,z,Rlow=S.symbols("C D G T B v z Rlow",commutative=False)
    conformal=-G*(T+B*v)+z
    original=C*conformal+D*v+Rlow
    reduced=(D-C*G*B)*v+C*z+Rlow-C*G*T
    zero("noncommutative_Schur_elimination",S.expand(original-reduced))


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
            conformal_coefficients()
            radial_and_indicial()
            coordinate_constraint_blocks()
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
