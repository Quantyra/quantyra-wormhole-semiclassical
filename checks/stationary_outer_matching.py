"""Focused exact checks supporting S139's fixed-reference outer theorem.

The ordinary analytic argument controls the boundary remainder. This script
checks algebra and normalizations; it does not certify a uniform heat bound.
Each execution is retained, including failed input snapshots and raw output.
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
RECEIPT = ROOT / "evidence/stationary-outer-matching.json"
LOG = ROOT / "evidence/stationary-outer-matching.log"
INPUTS = [
    "checks/stationary_outer_matching.py",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-exterior-boundary-weyl.md",
    "notes/stationary-massive-boundary-layer.md",
    "notes/stationary-paired-dirac-phase.md",
    "notes/spherical-anomaly-action.md",
    "notes/spherical-anomaly-sources.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-lll-quantum-exterior.md",
    "notes/coupled-fermion-domain.md",
]
CONTRACT = {
    "acceptance": "all exact symbolic residuals zero; all scope guards true",
    "numerical_tolerances": "none: no floating-point acceptance tests",
    "scope": "fixed RN reference, identical paired static domain, one-loop common scheme",
    "analytic_remainder": "ordinary scalar Brownian-bridge comparison; finite fixed-reference constants only",
    "excluded": ["uniform temperature constants", "full corrected exterior", "finite boundary charge",
                 "driven response transfer", "physical S074-S079 completion"],
    "failure_policy": "append every run and raw log; snapshot all failed-run input bytes",
}
ROWS = []


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    if not row["passed"]:
        print("FAIL " + json.dumps(row, sort_keys=True), flush=True)


def zero(label, value):
    value = S.factor(S.together(value))
    check(label, value == 0, value)


def bulk():
    r, b, lam, m, Q, q, G = S.symbols("r b lambda m Q q G", positive=True)
    f, a = S.Function("F")(r), S.Function("alpha")(r)
    sig = S.log(r/b)
    T = r*(S.diff(f,r)+f*S.diff(a,r))
    R = (-r*r*S.diff(f,r,2)-3*r*r*S.diff(a,r)*S.diff(f,r)
         -2*r*r*f*(S.diff(a,r,2)+S.diff(a,r)**2)+2*T-2*f)
    lag = S.exp(a)/r**2*(sig*(R**2/60-R/18+S.Rational(1,15)-q*q/3)
          +R*T/30-T*T/10+11*T*f/90-11*f*f/180-f/18)
    rn = 1+lam*r*r-2*m/r+Q/r**2
    def euler(z):
        return sum((-1)**j*S.diff(S.diff(lag,S.diff(z,r,j)),r,j) for j in range(3))
    ef = 2*(6*Q*sig-8*Q+3*m*r)/(15*r**4)
    ea = ((432*Q**2+288*Q*lam*r**4-864*Q*m*r+432*Q*r**2-60*q*q*r**4)*sig
          -593*Q**2-338*Q*lam*r**4+1408*Q*m*r-504*Q*r**2
          -33*lam**2*r**8+144*lam*m*r**5-588*m*m*r*r+288*m*r**3-11*r**4)/(180*r**6)
    zero("bulk_F_Euler", euler(f).subs({a:0,f:rn}).doit()-ef)
    zero("bulk_alpha_Euler", euler(a).subs({a:0,f:rn}).doit()-ea)
    rho, p = -ea/(4*S.pi), ea/(4*S.pi)-rn*ef/(2*S.pi)
    anomaly = (S.Rational(14,15)*(m/r**3-Q/r**4)**2+S.Rational(11,45)*Q**2/r**8
               -S.Rational(11,15)*lam*lam-q*q/(3*r**4))/(16*S.pi**2)
    jx = 2*G/r**2*(4*S.pi*r*r*anomaly-(ea-rn*ef)/(2*S.pi))
    zero("bulk_joint_trace", -rho+p+r*r*jx/(2*G)-4*S.pi*r*r*anomaly)
    zero("bulk_radial_Ward", S.diff(p,r)+S.diff(rn,r)*(rho+p)/(2*rn)-r*jx/(2*G))
    lapse = G/(15*S.pi)*(m/r**3+Q*(S.Rational(3,2)*sig-S.Rational(13,8))/r**4)
    zero("bulk_integrated_lapse", S.diff(lapse,r)+G*ef/(2*S.pi*r))
    xi = S.symbols("xi", positive=True)
    zero("bulk_extremal_lapse", lapse.subs({r:b,m:b*(1+2*xi),Q:b*b*(1+3*xi)})
         +G*(5+23*xi)/(120*S.pi*b*b))
    zero("bulk_vacuum_density", S.limit(rho/r**2,r,S.oo)-11*lam**2/(240*S.pi))
    zero("bulk_mass_log_coefficient", S.limit(r*(rho-11*lam**2*r*r/(240*S.pi)),r,S.oo)
         +lam*m/(5*S.pi))
    zero("bulk_lapse_tail", S.limit(r**3*lapse,r,S.oo)-G*m/(15*S.pi))
    zero("bulk_vacuum_pressure", S.limit(p/r**2,r,S.oo)+11*lam**2/(240*S.pi))
    print("bulk exact identities complete", flush=True)


def boundary():
    a, ap, app, s, sp, spp, Y, sig, t = S.symbols("a ap app s sp spp Y sigma t")
    def dx(expr):
        return S.diff(expr,a)*ap+S.diff(expr,ap)*app+S.diff(expr,s)*sp+S.diff(expr,sp)*spp+S.diff(expr,sig)*s
    r2 = -2*(ap+a*a)
    r0 = r2+2*Y
    rt = r0-6*t*(sp+a*s)-6*t*t*s*s
    kt, ks = a+t*s, t*s
    K, K2, K3 = kt+2*ks, kt**2+2*ks**2, kt**3+2*ks**3
    nt = -(ap+a*a)-t*(sp+a*s)
    ns, ts, ss = -t*sp, -t*a*s-t*t*s*s, Y-t*t*s*s
    curvature = (4*(-nt-2*ns)*K+18*(-nt*kt-2*ns*ks)
                 +4*(-2*ts*(kt+ks)-2*ss*ks))
    fpart = (3*(dx(rt)-2*t*s*rt)+5*rt*K+curvature
             +(34*K**3+78*K*K2-232*K3)/21)
    smear = -s*(36*K*K-54*K2)/7-6*K*(sp-t*s*s)
    P = [-6*app+S.Rational(44,7)*a**3+10*a*Y,
         -18*spp-18*a*sp-6*ap*s+12*a*a*s+10*Y*s,14*a*s*s,14*s**3]
    zero("boundary_full_f_polynomial", fpart-sum(P[j]*t**j for j in range(4)))
    Fcal = sig*(P[0]+P[1]/2+P[2]/3+P[3]/4)+(S.Rational(18,7)*a*a*s-6*a*sp
            -S.Rational(33,7)*a*s*s-9*s*sp-S.Rational(12,7)*s**3)
    zero("boundary_literal_path", S.integrate(sig*fpart+smear,(t,0,1))-Fcal)
    Jcal = s*(r0-3*sp-3*a*s-2*s*s)-sig*(dx(r0)-3*spp-3*ap*s-3*a*sp-s*r0+4*a*s*s+3*s**3)
    literalJ = S.integrate(rt*s-sig*(dx(rt)-2*t*s*rt),(t,0,1))
    zero("boundary_BoxR_conversion", literalJ-Jcal)
    A, d = S.Rational(11,360), S.Rational(7,540)
    Bflux = d*(r0*s-sig*dx(r2))-2*A*(sig*(spp+ap*s+a*sp)-(sp+a*s)*s+(r2-4*Y)*sig*s/3)
    Ccal = -Fcal-S.Rational(7,6)*Jcal+90*Bflux+3*(sig*dx(r2)-s*r2)
    Gcal = sig*(12*a*ap+S.Rational(44,7)*a**3+10*a*Y)-6*ap*s-S.Rational(24,7)*a*a*s-6*a*sp-S.Rational(96,7)*a*s*s-18*s*sp-S.Rational(85,21)*s**3
    zero("boundary_all_conversion_fluxes", Ccal+Gcal)
    k = S.symbols("k")
    zero("boundary_flat_ball_cubic", (34*(3*k)**3+78*(3*k)*(3*k*k)-232*(3*k**3))/21-44*k**3)
    zero("boundary_flat_ball_a4", 4*S.Rational(1,360)/(4*S.pi)**2*2*S.pi**2*44-S.Rational(11,180))
    u, X, am, hm, b = S.symbols("u X am hm b", positive=True)
    B = X+u*u-2*am*u**3+hm*u**4
    bp, bpp, bppp = (S.diff(B,u,j) for j in (1,2,3))
    jets = {a:-bp/(2*b*S.sqrt(B)),s:S.sqrt(B)/(b*u),Y:1/b**2,
            ap:(bpp/2-bp**2/(4*B))/b**2,
            app:-S.sqrt(B)*(bppp-bp*bpp/B+bp**3/(2*B**2))/(2*b**3),
            sp:(B/u**2-bp/(2*u))/b**2,
            spp:S.sqrt(B)*(2*B/u**3-3*bp/(2*u*u)+bpp/(2*u))/b**3}
    Gconst = S.factor(S.sqrt(B)*b**3*Gcal.subs(sig,0).subs(jets))
    Gidentity = -3*B*bpp/u-S.Rational(6,7)*bp*bp/u+S.Rational(132,7)*B*bp/u**2-S.Rational(463,21)*B*B/u**3
    zero("boundary_total_RN_exact", Gconst-Gidentity)
    zero("boundary_total_RN_finite", S.limit(Gconst+S.Rational(463,21)*X*X/u**3+S.Rational(260,21)*X/u,u,0)-S.Rational(232,21)*X*am)
    Fsig = S.factor(S.sqrt(B)*b**3*S.diff(Fcal,sig).subs(jets))
    Fconst = S.factor(S.sqrt(B)*b**3*Fcal.subs(sig,0).subs(jets))
    zero("boundary_literal_no_mass_log", S.limit(Fsig+S.Rational(29,2)*X*X/u**3+S.Rational(14,3)*X/u,u,0))
    zero("boundary_literal_RN_finite", S.limit(Fconst+S.Rational(75,7)*X*X/u**3+S.Rational(12,7)*X/u,u,0)+S.Rational(114,7)*X*am)
    zero("boundary_total_no_log", S.limit(S.sqrt(B)*b**3*S.diff(Gcal,sig).subs(jets),u,0))
    print("boundary exact identities complete", flush=True)


def paired():
    I = S.eye(2)
    p1, p2, p3 = S.Matrix([[0,1],[1,0]]), S.Matrix([[0,-S.I],[S.I,0]]), S.diag(1,-1)
    kron = S.kronecker_product
    gam = [kron(I,p1,I),-kron(p3,p3,I),-kron(p3,p2,p1),-kron(p3,p2,p2)]
    J, T3 = kron(p1,p1,I), kron(p3,I,I)
    Vgen = S.I*kron(p2,p1,I)
    chi = T3*gam[0]
    def mz(label,M):
        check(label,M==S.zeros(8),M if M!=S.zeros(8) else "zero matrix")
    for i in range(4):
        for j in range(i,4):
            mz(f"Clifford_{i}_{j}",gam[i]*gam[j]+gam[j]*gam[i]-2*int(i==j)*S.eye(8))
        mz(f"J_commutes_gamma_{i}",J*gam[i]-gam[i]*J)
        mz(f"T3_commutes_gamma_{i}",T3*gam[i]-gam[i]*T3)
        mz(f"constant_angle_commutes_gamma_{i}",Vgen*gam[i]-gam[i]*Vgen)
    mz("J_unitary",J.H*J-S.eye(8))
    mz("J_hermitian",J.H-J)
    mz("J_flips_chi",J*chi+chi*J)
    mz("T3_preserves_domain",T3*chi-chi*T3)
    mz("T3_anticommutes_J",T3*J+J*T3)
    mz("Hamiltonian_radial_orientation",gam[1]*gam[0]+S.I*kron(p3,p2,I))
    mz("Hamiltonian_angular_orientation",S.I*gam[1]*gam[2]-kron(I,p1,p1))
    plus, minus = (S.eye(8)+chi)/2,(S.eye(8)-chi)/2
    mz("J_maps_complementary_domains",J*plus-minus*J)
    check("balanced_projector_ranks",plus.trace()==minus.trace()==4)
    print("paired Clifford identities complete", flush=True)


def massive():
    x,y,v,z,b,c,n,q,r,G,lam,m = S.symbols("x y v z b c n q r G lambda m",positive=True)
    integ = S.integrate((x-y)**3/x**5,(x,y,S.oo))
    zero("cubic_Tonelli_integral",6*v*v*integ-3*v*v/(2*y))
    moment = S.integrate((2*z*(1-z))**2,(z,0,1))
    zero("Brownian_bridge_moment",moment-S.Rational(2,15))
    zero("cubic_angular_trace_coefficient",c*b*b/S.pi*3*moment/(2*y)-c*b*b/(5*S.pi*y))
    zero("angular_summable_tail_identity",2*(q+2*n)/(n*n*(n+q)**2)-2/q*(1/n**2-1/(n+q)**2))
    tau, F, I = S.Function("tau")(r),S.Function("F")(r),S.Function("I")(r)
    B = b*b*F/r**2
    p, rho = I/(2*F),I/(2*F)-b*b*tau/r**2
    jx = -2*G*b*b*tau/r**4
    ward = S.diff(p,r)+S.diff(F,r)*(rho+p)/(2*F)-r*jx/(2*G)
    zero("massive_reconstructed_Ward",ward.subs(S.diff(I,r),S.diff(B,r)*tau))
    tnull = -B*tau/4+I/4
    zero("massive_conservation_by_parts",(S.diff(tnull,r)+B*S.diff(tau,r)/4).subs(S.diff(I,r),S.diff(B,r)*tau))
    zero("massive_trace_dictionary",p-rho-b*b*tau/r**2)
    A = -m*lam/(5*S.pi*b*b)
    zero("mass_tail_cancellation",-b*b*A-m*lam/(5*S.pi))
    lapsehat = -G*m/(15*S.pi*r**3)
    zero("massive_lapse_leading_constraint",S.diff(lapsehat,r)-G*m/(5*S.pi*r**4))
    zero("lapse_tail_cancellation",lapsehat+G*m/(15*S.pi*r**3))
    t = S.symbols("t",positive=True)
    zero("log_remainder_small_time_split",S.integrate(t**-1*(S.sqrt(t)/y)**4,(t,0,y*y))-S.Rational(1,2))
    zero("integrable_outer_majorant",S.integrate((1+S.log(r))/r**2,(r,1,S.oo))-2)
    zero("AdS_radius_center_equation",(3*(lam-S.symbols("lambda0"))/(8*S.pi*G)+11*lam*lam/(960*S.pi**2))*8*S.pi*G/3
         -(lam-S.symbols("lambda0")+11*G*lam*lam/(360*S.pi)))
    print("massive boundary and reconstruction identities complete", flush=True)


def hashes(data):
    return {name:hashlib.sha256(value).hexdigest() for name,value in data.items()}


def main():
    before = {name:(ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
           "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
           "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    exit_code = 0
    with LOG.open("a",encoding="utf-8",newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN " + run["started_utc"],flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            for task in (bulk,boundary,paired,massive):
                task()
            check("scope_excludes_full_exterior","full corrected exterior" in CONTRACT["excluded"])
            check("scope_fixed_reference_constants","fixed-reference" in CONTRACT["analytic_remainder"])
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
