"""Exact algebra supporting the S139 solver and charge route correction.

The localized unboundedness proof and the full-action charge counterexample
are ordinary arguments. This run checks their curvature and normalization
identities; it does not compute the physical Dirac boundary charge.
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
RECEIPT = ROOT / "evidence/stationary-solver-charge.json"
LOG = ROOT / "evidence/stationary-solver-charge.log"
INPUTS = [
    "checks/stationary_solver_charge.py",
    "notes/stationary-exterior-solver-obstruction.md",
    "notes/stationary-exterior-charge-audit.md",
    "notes/spherical-anomaly-action.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-massive-boundary-layer.md",
]
CONTRACT = {
    "acceptance": "all symbolic residuals exactly zero; all scope guards true",
    "numerical_tolerances": "none",
    "scope": "retained one-Dirac anomaly scheme; local RN annulus; separate f(R) counterexample",
    "ordinary_arguments": "localized high-frequency norm contradiction; covariant relative charge scaling",
    "excluded": ["full stationary stress bound", "physical Dirac boundary charge",
                 "global Planck mode", "nonexistence of corrected exterior", "physical completion"],
    "failure_policy": "append every execution and raw log; retain failed input bytes",
}
ROWS = []


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    if not row["passed"]:
        print("FAIL " + json.dumps(row, sort_keys=True), flush=True)


def zero(label, expr):
    expr = S.factor(S.trigsimp(S.together(expr)))
    check(label, expr == 0, expr)


def curvature():
    U, r, th, ph = S.symbols("U r theta phi", real=True)
    coords = (U, r, th, ph)
    F, alpha, v = (S.Function(n)(r) for n in ("F", "alpha", "v"))
    diagonal = (-S.exp(2*alpha)*F, 1/F, r*r, r*r*S.sin(th)**2)
    metric = S.diag(*diagonal)
    inv = S.diag(*(1/x for x in diagonal))
    gamma = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                gamma[k,i,j] = S.simplify(sum(inv[k,l]*(
                    S.diff(metric[l,j],coords[i])+S.diff(metric[l,i],coords[j])
                    -S.diff(metric[i,j],coords[l])) for l in range(4))/2)
    ricci = []
    for i in range(4):
        value = sum(S.diff(gamma[k,i,i],coords[k])-S.diff(gamma[k,i,k],coords[i])
                    +sum(gamma[k,i,i]*gamma[l,k,l]-gamma[l,i,k]*gamma[k,i,l]
                         for l in range(4)) for k in range(4))
        ricci.append(S.simplify(value))
    direct_R = S.simplify(sum(inv[i,i]*ricci[i] for i in range(4)))
    expected_R = (-S.diff(F,r,2)-3*S.diff(alpha,r)*S.diff(F,r)
                  -2*F*(S.diff(alpha,r,2)+S.diff(alpha,r)**2)
                  -4*S.diff(F,r)/r-4*F*S.diff(alpha,r)/r+2*(1-F)/r**2)
    zero("four_dimensional_Christoffel_curvature",direct_R-expected_R)
    eps = S.symbols("epsilon")
    delta_R = S.diff(direct_R.subs(alpha,eps*v).doit(),eps).subs(eps,0)
    expected_delta = -2*F*S.diff(v,r,2)-(3*S.diff(F,r)+4*F/r)*S.diff(v,r)
    zero("lapse_curvature_variation",delta_R-expected_delta)
    lam, m, Q, q, G = S.symbols("lambda m Q q G", positive=True)
    rn = 1+lam*r*r-2*m/r+Q/r**2
    R0 = direct_R.subs({alpha:0,F:rn}).doit()
    zero("RN_constant_scalar_curvature",R0+12*lam)
    zero("RN_scalar_gauge_variation",S.diff(R0,r))
    u = S.Function("u")(r)
    density = S.exp(alpha)*r*r*S.sin(th)
    box = S.diff(density*F*S.diff(u,r),r)/density
    zero("radial_box_from_volume",box-(F*S.diff(u,r,2)
         +(S.diff(F,r)+F*S.diff(alpha,r)+2*F/r)*S.diff(u,r)))
    box0 = lambda z: F*S.diff(z,r,2)+(S.diff(F,r)+2*F/r)*S.diff(z,r)
    variation_box_R = box0(delta_R)
    zero("fourth_derivative_coefficient",S.expand(variation_box_R).coeff(S.diff(v,r,4))+2*F**2)
    zero("constant_background_delta_box_term",S.diff(R0,r,2)+S.diff(R0,r))
    magnetic_square = 2*(q*S.sin(th)/2)**2*inv[2,2]*inv[3,3]
    zero("magnetic_invariant",magnetic_square-q*q/(2*r**4))
    zero("magnetic_lapse_variation",S.diff(magnetic_square.subs(alpha,eps*v),eps))
    coefficient = 8*S.pi*G*S.Rational(1,30)/(16*S.pi**2)
    zero("trace_box_coefficient",coefficient-G/(60*S.pi))
    p, n, F0 = S.symbols("p n F0", positive=True)
    zero("proper_radial_symbol",(1-coefficient*F0*n*n).subs(n,p/S.sqrt(F0))
         -(1-G*p*p/(60*S.pi)))
    zero("local_principal_root",(1-G*p*p/(60*S.pi)).subs(p*p,60*S.pi/G))
    # The preconditioned trace update is -coefficient*Box, hence +coefficient*p^2.
    zero("Einstein_trace_update_multiplier",-coefficient*(-F0*n*n)-G*F0*n*n/(60*S.pi))
    cR, R, d, gradient2 = S.symbols("cR R box_sigma gradient_sigma_squared")
    weyl_density = cR*(R-6*eps*d-6*eps**2*gradient2)**2
    zero("finite_R2_Weyl_variation",S.diff(weyl_density,eps).subs(eps,0)+12*cR*R*d)
    zero("finite_R2_scheme_shift",16*S.pi**2*(-12*cR)+192*S.pi**2*cR)
    print("curvature and retained-scheme identities complete",flush=True)


def charge():
    R, Lam, beta, G, L, M, x = S.symbols("R Lambda beta G L M x", positive=True)
    f = R-2*Lam+beta*R**2
    fprime = S.diff(f,R)
    einstein_coefficient = fprime*R/4-f/2
    zero("fR_constant_curvature_equation",einstein_coefficient-(Lam-R/4))
    zero("R2_Einstein_bulk_variation",2*R*(R/4)-R**2/2)
    zero("fR_Einstein_family",einstein_coefficient.subs(R,4*Lam))
    phi0 = fprime.subs(R,4*Lam)
    zero("fR_effective_Einstein_factor",phi0-(1+8*beta*Lam))
    zero("fR_AdS_charge_factor",fprime.subs(R,-12/L**2)-(1-24*beta/L**2))
    cR = S.symbols("cR")
    zero("fR_coupling_conversion",fprime.subs(beta,16*S.pi*G*cR)-(1+32*S.pi*G*cR*R))
    zero("constant_curvature_action_scaling",(f-phi0*(R-2*Lam)).subs(R,4*Lam))
    # Direct series of the stated Einstein Brown-York expression at r=1/x.
    sqrtF = S.sqrt(1+L**2*x**2-2*G*M*L**2*x**3)/(L*x)
    by = sqrtF/G*(1/(L*x*x)+L/2-sqrtF/x)
    zero("Einstein_Brown_York_finite_mass",S.series(by,x,0,1).removeO()-M)
    offset = S.symbols("source_only_offset")
    M0 = S.symbols("M0")
    zero("relative_charge_source_offset_cancellation",
         ((phi0*M+offset)-(phi0*M0+offset))-phi0*(M-M0))
    # Canonical-frame lapse/F vertices, componentwise on radial test spinor data.
    aa, FF, FFp, j, deltaF, deltaFp = S.symbols("a F Fprime j deltaF deltaFprime",positive=True)
    psip, psit, psi = S.symbols("psi_r psi_U psi")
    eps, eta = S.symbols("epsilon eta")
    radial = S.sqrt(FF)*psip+FFp*psi/(4*S.sqrt(FF))
    time = psit/(aa*S.sqrt(FF))
    zero("canonical_lapse_vertex",S.diff(time.subs(aa,aa*S.exp(eps*j)),eps).subs(eps,0)+j*time)
    zero("canonical_radial_F_vertex",S.diff(radial.subs({FF:FF+eps*deltaF,FFp:FFp+eps*deltaFp}),eps).subs(eps,0)
         -(deltaF*psip/(2*S.sqrt(FF))+deltaFp*psi/(4*S.sqrt(FF))-FFp*deltaF*psi/(8*FF**S.Rational(3,2))))
    varied_time = time.subs({aa:aa*S.exp(eta*j),FF:FF+eps*deltaF})
    zero("canonical_mixed_seagull",S.diff(varied_time,eps,eta).subs({eps:0,eta:0})
         -j*deltaF*psit/(2*aa*FF**S.Rational(3,2)))
    print("charge comparison and regulated vertex identities complete",flush=True)


def hashes(data):
    return {name:hashlib.sha256(value).hexdigest() for name,value in data.items()}


def main():
    before = {name:(ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
           "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
           "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+run["started_utc"],flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            curvature()
            charge()
            check("scope_excludes_physical_charge","physical Dirac boundary charge" in CONTRACT["excluded"])
            check("scope_excludes_global_mode","global Planck mode" in CONTRACT["excluded"])
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"] = hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        exit_code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if exit_code else "passed"
        run["terminal_exit_status"] = exit_code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if exit_code:
            for name,data in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(data).hexdigest(),"base64":base64.b64encode(data).decode("ascii")},sort_keys=True))
        print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),"terminal_exit_status":exit_code},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),"receipt":str(RECEIPT)},sort_keys=True))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
