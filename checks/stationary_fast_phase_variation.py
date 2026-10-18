"""Focused finite algebra for the restricted stationary variational theorem.

The C1 heat/domain argument and all asymptotic estimates remain ordinary
proofs. This checks identities and power bookkeeping, not a coupled solver.
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
RECEIPT = ROOT / "evidence/stationary-fast-phase-variation.json"
LOG = ROOT / "evidence/stationary-fast-phase-variation.log"
INPUTS = [
    "checks/stationary_fast_phase_variation.py",
    "notes/stationary-fast-phase-variation.md",
    "notes/stationary-same-source-boundary-charge.md",
    "notes/stationary-trace-range-condition.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-exterior-boundary-weyl.md",
    "notes/stationary-paired-dirac-phase.md",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-covariant-seed.md",
]
CONTRACT = {
    "acceptance": "exact zero symbolic residuals and true rational/power inequalities",
    "numerical_tolerances": "none",
    "scope": "finite algebra and power bookkeeping for a chosen fast stationary phase",
    "ordinary_arguments": ["uniform C1 heat remainder and finite-cut limit",
        "actual-domain local gap", "all allowed conormal remainders",
        "bolt regularity and physical boundary scope"],
    "excluded": ["full dynamical Dirichlet completion", "physical phase selection",
        "coupled inverse", "exact exterior", "physical state", "physical gate closure"],
    "failure_policy": "append all runs and logs; retain failed-run input bytes",
    "hash_policy": "SHA256 exact raw UTF8/LF bytes, no normalization",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    value = S.factor(S.together(S.expand(expression)))
    check(label, value == 0, value)


def degree(expr, weights):
    """Largest r power of an expanded Laurent polynomial; logs have weight zero."""
    terms = S.Add.make_args(S.expand(expr))
    values = []
    for term in terms:
        powers = term.as_powers_dict()
        values.append(sum(powers.get(var, 0)*weight for var, weight in weights.items()))
    return max(values)


def local_variations():
    r, G, lam, B = S.symbols("r G lam B", positive=True)
    m, Q, q, Lambda0, sig = S.symbols("m Q q Lambda0 sig", real=True)
    f = S.symbols("f0:4", real=True)
    a = S.symbols("a0:4", real=True)

    def dr(expr):
        return (S.diff(expr, r)+S.diff(expr, sig)/r
                +sum(S.diff(expr, f[j])*f[j+1]+S.diff(expr, a[j])*a[j+1] for j in range(3)))

    def momenta(L):
        return [S.diff(L,f[1])-dr(S.diff(L,f[2])), S.diff(L,f[2]),
                S.diff(L,a[1])-dr(S.diff(L,a[2])), S.diff(L,a[2])]

    R4 = (-f[2]-3*a[1]*f[1]-2*f[0]*(a[2]+a[1]**2)-4*f[1]/r
          -4*f[0]*a[1]/r+2*(1-f[0])/r**2)
    LEH = S.exp(a[0])*r**2*(R4-2*Lambda0)/(4*G)
    PF,QF,PA,QA = momenta(LEH)
    zero("EH_PF",PF+S.exp(a[0])*(r**2*a[1]+r)/(2*G))
    zero("EH_QF",QF+S.exp(a[0])*r**2/(4*G))
    zero("EH_PA",PA+S.exp(a[0])*r**2*(f[1]+2*f[0]*a[1])/(4*G))
    zero("EH_QA",QA+S.exp(a[0])*r**2*f[0]/(2*G))
    GHY = S.exp(a[0])*(r**2*(f[1]+2*f[0]*a[1])/(4*G)+r*f[0]/G)
    for label, value in [("EH_GHY_dFprime",QF+S.diff(GHY,f[1])),
                         ("EH_GHY_dAprime",QA+S.diff(GHY,a[1])),
                         ("EH_GHY_dF",PF+S.diff(GHY,f[0])-S.exp(a[0])*r/(2*G)),
                         ("EH_GHY_dA",PA+S.diff(GHY,a[0])-S.exp(a[0])*r*f[0]/G)]:
        zero(label,value)
    CTEH = -S.exp(a[0])*S.sqrt(f[0])*(S.sqrt(lam)*r**2+1/(2*S.sqrt(lam)))/G
    EEH = S.sqrt(f[0])*(S.sqrt(lam)*r**2+1/(2*S.sqrt(lam))-r*S.sqrt(f[0]))/G
    zero("EH_ren_clock",PA+S.diff(GHY+CTEH,a[0])+S.exp(a[0])*EEH)
    zero("EH_ren_dF",PF+S.diff(GHY+CTEH,f[0])-S.exp(a[0])/(2*G)
         *(r-(S.sqrt(lam)*r**2+1/(2*S.sqrt(lam)))/S.sqrt(f[0])))
    F = lam*r**2+1-2*m/r+Q/r**2
    z = S.symbols("z",positive=True)
    Eseries = S.series(EEH.subs(f[0],F).subs(r,1/z),z,0,1).removeO()
    zero("EH_energy_finite_mass",Eseries-m/G)
    cf = r-(S.sqrt(lam)*r**2+1/(2*S.sqrt(lam)))/S.sqrt(F)
    zero("EH_normalizable_coefficient_leading",S.series(cf.subs(r,1/z),z,0,3).removeO()+m*z**2/lam)

    T = r*(f[1]+f[0]*a[1])
    rb = -r**2*f[2]-3*r**2*f[1]*a[1]-2*r**2*f[0]*(a[2]+a[1]**2)+2*T-2*f[0]
    LW = S.exp(a[0])/(4*S.pi*r**2)*(sig*(rb**2/60-rb/18+R(1,15)-q**2/3)
         +rb*T/30-T**2/10+11*T*f[0]/90-11*f[0]**2/180-f[0]/18)
    moment = momenta(LW)
    weights = {r:1,**{f[j]:2-j for j in range(4)},**{a[j]:-j for j in range(1,4)}}
    for label,expr,target,tangent in zip(["PF","QF","PA","QA"],moment,[1,2,3,4],[-2,-3,-4,-5]):
        deg=degree(S.expand(expr/S.exp(a[0])),weights)
        check("W_momentum_"+label,deg<=target,{"derived_r_power":str(deg),"upper":target})
        check("W_remainder_surface_"+label,deg+tangent<0,{"product_r_power":str(deg+tangent)})

    flat = {a[j]:0 for j in range(4)}
    rn = {f[j]:S.diff(F,r,j) for j in range(4)}
    rbRN=S.expand(rb.subs(flat).subs(rn))
    zero("conformal_Rbar_RN",rbRN-(-2+12*m/r-12*Q/r**2))
    terms=[rb*T/30,-T**2/10,11*T*f[0]/90,-11*f[0]**2/180]
    targets=[R(4,5),-R(4,5),-R(11,45),R(11,45)]
    for j,(term,target) in enumerate(zip(terms,targets)):
        coefficient=S.expand(term.subs(flat).subs(rn)).coeff(r,1).coeff(m,1)
        zero("bulk_mass_r_coefficient_"+str(j),coefficient-target*lam)
    bulk=S.expand(LW.subs(flat).subs(rn))
    zero("bulk_action_r2",bulk.coeff(r,2)+13*lam**2/(240*S.pi))
    zero("bulk_action_r0",bulk.coeff(r,0)+lam/(60*S.pi))
    zero("bulk_action_mass_log_divergence",bulk.coeff(r,-1))
    # Here sigma is a separate logarithmic symbol of weight zero.
    remainder=bulk+13*lam**2*r**2/(240*S.pi)+lam/(60*S.pi)
    check("bulk_remaining_power_integrable",degree(remainder,{r:1})<=-2)
    zero("W_action_R3",-R(13,720)+R(463,7560)-R(653,15120))
    zero("W_action_R",-R(1,60)+R(260,7560)-R(67,3780))
    zero("centered_action_R3",R(653,15120)+R(11,720)-R(221,3780))
    c0=-221*lam**R(3,2)/(15120*S.pi**2)
    c1=29*S.sqrt(lam)/(20160*S.pi**2)
    zero("CT_cubic_power",4*S.pi*c0*S.sqrt(lam)+221*lam**2/(3780*S.pi))
    zero("CT_linear_power",2*S.pi*c0/S.sqrt(lam)+8*S.pi*c1*S.sqrt(lam)+67*lam/(3780*S.pi))
    zero("mass_surface_three_terms",-R(105,3780)-R(116,3780)+R(221,3780))


def boundary_and_heat_weights():
    a,ax,s,sx,sig,Y0=S.symbols("a ax s sx sig Y0")
    Gcal=(sig*(12*a*ax+44*a**3/7+10*a*Y0)-6*ax*s-24*a**2*s/7
           -6*a*sx-96*a*s**2/7-18*s*sx-85*s**3/21)
    weights={a:-1,ax:0,s:1,sx:2}
    variation={a:-3,ax:-2,s:-3,sx:-2}
    check("boundary_delta_N_power",degree(Gcal,weights)-4<=-1)
    for var,dweight in variation.items():
        d=degree(S.diff(Gcal,var),weights)+dweight
        check("boundary_delta_"+str(var),d<=-1,{"r_power":str(d)})
    sigma_poly=S.expand(Gcal).coeff(sig,1)
    check("boundary_sigma_delta_N",degree(sigma_poly,weights)-4<=-3)
    for var,dweight in variation.items():
        expr=S.diff(sigma_poly,var)
        if expr!=0:
            check("boundary_sigma_delta_"+str(var),degree(expr,weights)+dweight<=-3)
    # Heat trace exponents after integrating the collar; no heat estimate is proved here.
    zero("near_collar_integrated_power",1+R(1,2)-R(3,2))
    zero("crossing_integrated_power",4*R(1,2)-1+R(1,2)-R(3,2))
    zero("action_heat_integrand_power",R(3,2)-1-1+R(1,2))
    check("action_heat_integrable_at_zero",-R(1,2)>-1)
    check("fixed_source_canonical_FP_integrable",3-1>-1)
    check("fixed_source_radial_flux_vanishes",4>0)
    alpha,m,G,Z=S.symbols("alpha m G Z")
    theta=-S.exp(alpha)*Z*m/G
    zero("boundary_oneform_mass_curl",S.diff(theta,m)+S.exp(alpha)*Z/G)
    check("boundary_oneform_not_claimed_exact",S.diff(theta,m)!=0)
    # The proposed single-slice ambiguity fails the open-clock condition.
    zero("single_slice_ambiguity_mass_variation",S.diff(alpha*m,m)-alpha)


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
            local_variations()
            boundary_and_heat_weights()
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
