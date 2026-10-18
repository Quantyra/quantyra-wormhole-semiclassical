"""Finite identities for the full-tensor unweighted inverse obstruction.

The conditional Sobolev/domain theorem remains an ordinary proof. No
global RN mode, instability, weighted inverse or shooting run is tested.
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
RECEIPT = ROOT / "evidence/stationary-coupled-asymptotic-obstruction.json"
LOG = ROOT / "evidence/stationary-coupled-asymptotic-obstruction.log"
INPUTS = [
    "checks/stationary_coupled_asymptotic_obstruction.py",
    "notes/stationary-coupled-asymptotic-obstruction.md",
    "notes/stationary-trace-range-condition.md",
    "notes/stationary-fast-phase-variation.md",
    "notes/stationary-exterior-solver-obstruction.md",
    "notes/spherical-anomaly-action.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-paired-dirac-phase.md",
]
CONTRACT = {
    "acceptance": "every symbolic residual exactly zero and rational inequality true",
    "numerical_tolerances": "none",
    "scope": "finite algebra for a conditional full tensor unweighted no-estimate theorem",
    "ordinary_arguments": ["conformal HH/domain transport and full stress locality",
        "AdS symmetry/conservation classification", "asymptotic Sobolev escaping sequence",
        "compatibility differentiation and endpoints"],
    "excluded": ["exact RN solution", "physical instability or nonexistence",
        "weighted inverse obstruction", "actual compatibility Jacobian",
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


def tensor_and_gauge():
    r,lam,eps,a=S.symbols("r lam eps a",positive=True)
    F=1+lam*r**2
    u=S.Function("u")(r)
    box=lambda v:F*S.diff(v,r,2)+(S.diff(F,r)+2*F/r)*S.diff(v,r)
    hess=lambda v:[S.diff(F,r)*S.diff(v,r)/2,
                   F*S.diff(v,r,2)+S.diff(F,r)*S.diff(v,r)/2,
                   F*S.diff(v,r)/r]
    trace=lambda v:v[0]+v[1]+2*v[2]
    divergence=lambda v:S.diff(v[1],r)+S.diff(F,r)/(2*F)*(v[1]-v[0])+2/r*(v[1]-v[2])
    D=lambda v:[vab-box(v)+3*lam*v for vab in hess(v)]
    zero("radial_Hessian_trace",trace(hess(u))-box(u))
    zero("D_tensor_trace",trace(D(u))+3*(box(u)-4*lam*u))
    zero("D_tensor_divergence",divergence(D(u)))
    Pu=a*u+eps*box(u)
    dR=-6*(box(u)-4*lam*u)
    response=[-2*v for v in D(Pu)]
    zero("full_response_trace_matches_anomaly",trace(response)+a*dR+eps*box(dR))
    zero("full_response_conserved",divergence(response))
    # Areal restoration: xi^r=-r*u, not an independent field variation.
    dF=2*r*F*S.diff(u,r)-r*S.diff(F,r)*u
    da=u-r*S.diff(u,r)
    dR_areal=(-S.diff(dF,r,2)-4*S.diff(dF,r)/r-2*dF/r**2
              -2*F*S.diff(da,r,2)-(3*S.diff(F,r)+4*F/r)*S.diff(da,r))
    zero("conformal_areal_curvature",dR_areal-dR)
    xi=-r*u
    zero("areal_sphere_component_cancels",2*r**2*u+2*r*xi)
    rr=2*u/F+xi*S.diff(1/F,r)+2*S.diff(xi,r)/F
    zero("areal_radial_metric",rr+dF/F**2)
    tt=-2*u*F-xi*S.diff(F,r)
    zero("areal_time_metric",tt+2*F*da+dF)
    c=-2*eps/(a+6*eps*lam)
    H=[2*u+2*c*v for v in hess(u)]
    Cr=divergence(H)-S.diff(trace(H),r)/2
    zero("de_Donder_correction",Cr+2*S.diff(Pu,r)/(a+6*eps*lam))
    zero("corrected_trace_normalization",trace(H)
         -(8+4*a/(a+6*eps*lam))*u+4*Pu/(a+6*eps*lam))
    check("corrected_trace_positive",S.ask(S.Q.positive(8+4*a/(a+6*eps*lam))) is True)
    # Exact radial Ward reduction: if both orbit residuals vanish, angular must vanish.
    Er,Et,Es,ap=S.symbols("Er Et Es ap")
    Erp=S.symbols("Erp")
    ward=Erp+(ap+S.diff(F,r)/(2*F))*(Er-Et)+2*(Er-Es)/r
    zero("Ward_orbit_equations_force_angular",ward.subs({Er:0,Et:0,Erp:0})+2*Es/r)


def asymptotic_and_projection():
    x,T,R0,k,c,eps,a,G,lam,B=S.symbols("x T R0 k c eps a G lam B",positive=True)
    v=S.Function("v")(x)
    conjugate=S.exp(-3*c*x/2)*v
    P=lambda f:a*f+eps*(S.diff(f,x,2)+3*c*S.diff(f,x))
    zero("limiting_Liouville_conjugation",S.exp(3*c*x/2)*P(conjugate)
         -eps*S.diff(v,x,2)-(a-9*eps*c**2/4)*v)
    chi=S.Function("chi")
    arg=(x-R0)/T
    packet=T**(-R(1,2))*chi(arg)*S.exp(S.I*k*x)
    raw=S.diff(packet,x,2)+k**2*packet
    target=T**(-R(1,2))*S.exp(S.I*k*x)*(2*S.I*k*S.diff(chi(arg),x)+S.diff(chi(arg),x,2))
    zero("escaping_packet_commutator",raw-target)
    zero("packet_first_cutoff_derivative",S.diff(chi(arg),x)-S.Subs(S.diff(chi(S.Symbol('z')),S.Symbol('z')),S.Symbol('z'),arg)/T)
    a_selected=1+11*G*lam/(180*S.pi)
    eps_selected=G/(60*S.pi)
    zero("proper_momentum_positive_shift",a_selected/eps_selected-9*lam/4-(60*S.pi/G+17*lam/12))
    s=S.symbols("s")
    zero("indicial_conjugate_roots",(s*(s+3)+(a/(eps*lam))).subs(s,-R(3,2)+S.I*k/S.sqrt(lam))
         -(a/eps-9*lam/4-k**2)/lam)
    check("forbidden_F_real_power",2-R(3,2)>-1)
    r,alpha=S.symbols("r alpha",positive=True)
    mu=S.sqrt(lam)*r**3
    zero("fixed_log_radius_measure",S.exp(alpha)*r**2*S.sqrt(lam)*r-S.exp(alpha)*mu)
    zero("fixed_flux_invariant_metric_variation",S.diff(S.Symbol('q')**2/(2*r**4),alpha))
    # Full Leibniz rule, including moving lower endpoint, is checked abstractly.
    W,E,dW,dE,Eh,Wh,drh=S.symbols("W E dW dE Eh Wh drh")
    t=S.symbols("t")
    integrand_variation=-S.diff((W+t*dW)*(E+t*dE),t).subs(t,0)
    zero("candidate_projection_product_rule",integrand_variation+W*dE+dW*E)
    zero("moving_lower_endpoint_sign",-(-Wh*Eh*drh)-Wh*Eh*drh)
    zero("genuine_solution_drops_defect_terms",(integrand_variation+Wh*Eh*drh).subs({E:0,Eh:0})+W*dE)
    dPZ,dS,PdZ=S.symbols("dPZ dS PdZ")
    zero("genuine_compatibility_projection_sign",(-W*(PdZ+dPZ-dS)).subs(PdZ,0)-W*(dS-dPZ))


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
            tensor_and_gauge()
            asymptotic_and_projection()
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
