"""Exact algebra for the conditional same-source boundary charge.

Distributional continuation and actual-domain limits remain reviewed
ordinary arguments. A full gravitational Dirichlet completion remains an
unproved hypothesis; this checker does not prove or physically select it.
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
RECEIPT = ROOT / "evidence/stationary-same-source-boundary-charge.json"
LOG = ROOT / "evidence/stationary-same-source-boundary-charge.log"
INPUTS = [
    "checks/stationary_same_source_boundary_charge.py",
    "notes/stationary-same-source-boundary-charge.md",
    "notes/stationary-exterior-charge-audit.md",
    "notes/stationary-exterior-boundary-weyl.md",
    "notes/stationary-exterior-conformal-domain.md",
    "notes/stationary-paired-dirac-phase.md",
    "notes/stationary-massive-boundary-layer.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-centered-conormal-source.md",
    "notes/uniform-dirac-reference.md",
]
CONTRACT = {
    "acceptance": "all symbolic residuals exactly zero and rational inequalities true",
    "numerical_tolerances": "none",
    "scope": "focused algebra of a charge conditional on a full open-source Dirichlet completion",
    "ordinary_arguments": [
        "meromorphic distribution and finite-cut fixed-test limit",
        "actual paired Dirac domain, phase, and no extra surface spectrum",
        "locally integrable remainder and total bulk cancellation before collar limit",
        "full Dirichlet completion hypothesis and source-only charge ambiguity",
    ],
    "excluded": [
        "formal analytic certificate", "physical selection of a gravitational boundary ensemble",
        "all-normal-jet completion", "exact exterior existence", "physical state matching",
        "pole normalization of Newton coupling", "physical gate closure",
    ],
    "failure_policy": "append every run/log and retain failed-run input bytes",
    "hash_policy": "SHA256 of exact raw bytes; UTF8/LF/noBOM; no normalization",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    value = S.simplify(S.factor(S.together(expression)))
    check(label, value == 0, value)


def mellin_and_distributions():
    s = S.symbols("s", real=True)
    v, y, t, mu, ell = S.symbols("v y t mu ell", positive=True)
    # Direct Gaussian: D=2 integral_y^infinity (x-y)^3 p_(vt)(x) dx.
    # The t integral is Gamma(5/2-s)*(x^2/(2v))^(s-5/2).
    gaussian_prefactor = 2/S.sqrt(2*S.pi*v)*(2*v)**(R(5,2)-s)
    zero("Gaussian_variance_two_prefactor", gaussian_prefactor
         -2**(3-s)*v**(2-s)/S.sqrt(S.pi))
    beta_x = S.gamma(4)*S.gamma(1-2*s)/S.gamma(5-2*s)
    raw_I = gaussian_prefactor*S.gamma(R(5,2)-s)*beta_x
    # Explicit duplication identity keeps the variance factors auditable.
    duplication = 2**(4-2*s)*S.gamma(R(5,2)-s)*S.gamma(3-s)/S.sqrt(S.pi)
    I_after_duplication = (gaussian_prefactor*S.gamma(R(5,2)-s)
                           *S.gamma(4)*S.gamma(1-2*s)/duplication)
    I_target = 3*2**s*v**(2-s)*S.gamma(1-2*s)/S.gamma(3-s)
    zero("Gaussian_Mellin_spatial_beta_after_duplication", I_after_duplication/I_target-1)
    check("duplication_identity_at_zero", S.simplify(duplication.subs(s,0)-S.gamma(5))==0)
    zero("raw_Mellin_normalization_at_zero", raw_I.subs(s,0)-R(3,2)*v**2)
    zero("raw_Mellin_log_derivative_at_zero",
         (S.diff(raw_I,s)/raw_I).subs(s,0)
         -(S.log(2)-S.log(v)+S.EulerGamma+R(3,2)))
    z = S.symbols("z", real=True)
    angular_v_beta = 2**(2-s)*S.gamma(3-s)**2/S.gamma(6-2*s)
    F_derived = S.simplify(I_target/v**(2-s)*angular_v_beta/S.pi)
    F = 12*S.gamma(1-2*s)*S.gamma(3-s)/(S.pi*S.gamma(6-2*s))
    zero("variance_factor_cancels_no_two_power", F_derived/F-1)
    zero("F_zero", F.subs(s,0)-1/(5*S.pi))
    zero("F_log_derivative", (S.diff(F,s)/F).subs(s,0)
         -(S.EulerGamma+R(46,15)))
    zero("beta_variance_fourth", S.integrate((2*z*(1-z))**2,(z,0,1))-R(2,15))
    gaussian_x = S.symbols("gaussian_x", positive=True)
    half_fourth = S.integrate(gaussian_x**4*S.exp(-gaussian_x**2/(2*v*t))
                              /S.sqrt(2*S.pi*v*t),(gaussian_x,0,S.oo))/2
    zero("Gaussian_half_fourth_moment", half_fourth-R(3,4)*(v*t)**2)
    c, B = S.symbols("c B", real=True)
    layer_integrated_y = c/(2*S.pi*t)*R(3,4)*t**2*R(2,15)
    zero("cubic_heat_profile_integrated_y", layer_integrated_y-c*t/(20*S.pi))
    a4_layer = S.simplify((2*B**2/t)*layer_integrated_y)
    A0 = c*B**2/(5*S.pi)
    zero("boundary_heat_residue", a4_layer-A0/2)
    normalized = mu**(2*s)*F/S.gamma(1+s)/F.subs(s,0)
    zero("mass_Green_meromorphic_surface_constant", S.diff(normalized,s).subs(s,0)/2
         -(S.log(mu)+S.EulerGamma+R(23,15)))
    # y^(2s-1): residue delta/(2s). A fixed FP length adds log(ell).
    Ccan = S.log(mu*ell)+S.EulerGamma+R(23,15)-R(1,2)
    zero("negative_boundary_anomaly_shift", Ccan
         -(S.log(mu*ell)+S.EulerGamma+R(31,30)))
    lam, r = S.symbols("lam r", positive=True)
    y_match = B/(S.sqrt(lam)*r)
    # Use positive B in logarithm identities below.
    Bp = S.symbols("Bp", positive=True)
    y_match = Bp/(S.sqrt(lam)*r)
    matching = S.expand_log(S.log(y_match/ell)+S.log(mu*ell),force=True)
    zero("matching_length_cancels", matching-(-S.log(r/Bp)+S.log(mu/S.sqrt(lam))))
    zero("local_constant_combination", R(31,30)-R(7,36)-R(151,180))
    # The local and canonical limiting distributions have opposite FP kernels.
    A, FP = S.symbols("A FP")
    Cw = -A*S.log(S.sqrt(lam)*ell)-A*R(7,36)
    total = A*FP+A*Ccan-A*FP+Cw
    zero("full_surface_distribution_and_scale", S.expand_log(total,force=True)
         -A*(S.log(mu/S.sqrt(lam))+S.EulerGamma+R(151,180)))
    # Proper lapse = orbit Weyl plus radial diffeomorphism.
    aa = S.Function("a")(y)
    jj = S.Function("j")(y)
    integral = S.Integral(jj/aa,(y,0,y))
    # Use an independent integration symbol to avoid a bound/free ambiguity.
    q = S.symbols("q", positive=True)
    integral = S.Integral(S.Function("j")(q)/S.Function("a")(q),(q,0,y))
    xi = -aa*integral
    omega = -S.diff(xi,y)
    zero("lapse_Weyl_radial_decomposition", omega+S.diff(aa,y)*xi/aa-jj)


def local_variation():
    r, B, lam, m, Q, G = S.symbols("r B lam m Q G", positive=True)
    f0,f1,f2,f3 = S.symbols("f0 f1 f2 f3")
    a0,a1,a2,a3 = S.symbols("a0 a1 a2 a3")
    sig = S.log(r/B)
    jet_map = ((f0,f1),(f1,f2),(f2,f3),(a0,a1),(a1,a2),(a2,a3))

    def dr(expr):
        return S.diff(expr,r)+sum(S.diff(expr,left)*right for left,right in jet_map)

    Rbar = -r**2*f2-3*r**2*f1*a1-2*r**2*f0*(a2+a1**2)+2*r*(f1+f0*a1)-2*f0
    T = r*(f1+f0*a1)
    qflux = S.symbols("qflux")
    LW = S.exp(a0)/(4*S.pi*r**2)*(sig*(Rbar**2/60-Rbar/18+R(1,15)-qflux**2/3)
         +Rbar*T/30-T**2/10+11*T*f0/90-11*f0**2/180-f0/18)
    flat_alpha = {a0:0,a1:0,a2:0,a3:0}
    Ja = (sig*(Rbar/30-R(1,18))+T/30).subs(flat_alpha)
    Ka = (Rbar/30-T/5+11*f0/90).subs(flat_alpha)
    Qa = S.diff(LW,a2).subs(flat_alpha)
    Pa = (S.diff(LW,a1)-dr(S.diff(LW,a2))).subs(flat_alpha)
    zero("bulk_alpha_second_momentum", Qa+f0*Ja/(2*S.pi))
    zero("bulk_alpha_first_momentum", Pa-((-f1+2*f0/r)*Ja+2*f0*dr(Ja)+f0*Ka/r)/(4*S.pi))
    Qf = S.diff(LW,f2).subs(flat_alpha)
    Pf = (S.diff(LW,f1)-dr(S.diff(LW,f2))).subs(flat_alpha)
    zero("bulk_F_second_momentum", Qf+Ja/(4*S.pi))
    zero("bulk_F_first_momentum", Pf-((2*Ja+Ka)/r+dr(Ja))/(4*S.pi))
    F = lam*r**2+1-2*m/r+Q/r**2
    rn = {f0:F,f1:S.diff(F,r),f2:S.diff(F,r,2),f3:S.diff(F,r,3)}
    PaRN = S.simplify(Pa.subs(rn))
    PfRN = S.simplify(Pf.subs(rn))
    QfRN = S.simplify(Qf.subs(rn))
    z = S.symbols("z", positive=True)

    def outer(expr):
        prepared = S.expand_log(expr,force=True).expand()
        prepared = S.expand_log(prepared.subs(r,1/z),force=True).expand()
        return S.series(prepared,z,0,1).removeO().expand()

    expected_Pa = -lam**2/(360*S.pi*z**3)-lam/(60*S.pi*z)
    expected_Pa += lam*m*S.log(B*z)/(5*S.pi)+lam*m/(90*S.pi)
    # The claimed expansion is through its finite part plus o(1).
    # SymPy retains z*log(z) at order 1 because it is not O(z).
    zero("bulk_clock_outer_series", S.limit(outer(PaRN)-expected_Pa,z,0,dir="+"))

    # Independent proper derivatives of the complete boundary polynomial.
    N = B*S.sqrt(F)/r
    dx = lambda expr: r*S.sqrt(F)*S.diff(expr,r)/B
    a = S.simplify(dx(N)/N)
    ss = S.sqrt(F)/B
    ax, sx = dx(a), dx(ss)
    Gcal = sig*(12*a*ax+44*a**3/7+10*a/B**2)-6*ax*ss-24*a**2*ss/7
    Gcal += -6*a*sx-96*a*ss**2/7-18*ss*sx-85*ss**3/21
    BW = S.simplify(-B**2*N*Gcal/(360*S.pi))
    expected_BW = 463*lam**2/(7560*S.pi*z**3)+260*lam/(7560*S.pi*z)
    expected_BW += -232*lam*m/(7560*S.pi)
    zero("full_boundary_Weyl_polynomial_series", outer(BW)-expected_BW)
    clock = S.simplify(PaRN+BW)
    clock_target = 221*lam**2/(3780*S.pi*z**3)+67*lam/(3780*S.pi*z)
    clock_target += lam*m*S.log(B*z)/(5*S.pi)-37*lam*m/(1890*S.pi)
    zero("full_local_clock_series", outer(clock)-clock_target)
    c0,c1 = S.symbols("c0 c1")
    Bct = 4*S.pi*S.sqrt(F)*(c0*r**2+2*c1)
    ct_series = outer(Bct)
    ct_target = 4*S.pi*c0*S.sqrt(lam)/z**3
    ct_target += (2*S.pi*c0/S.sqrt(lam)+8*S.pi*c1*S.sqrt(lam))/z
    ct_target += -4*S.pi*c0*m/S.sqrt(lam)
    zero("intrinsic_counterterm_covariant_series", ct_series-ct_target)
    chosen = S.solve([4*S.pi*c0*S.sqrt(lam)+221*lam**2/(3780*S.pi),
                     2*S.pi*c0/S.sqrt(lam)+8*S.pi*c1*S.sqrt(lam)+67*lam/(3780*S.pi)],(c0,c1))
    zero("required_volume_counterterm", chosen[c0]+221*lam**R(3,2)/(15120*S.pi**2))
    zero("required_curvature_counterterm", chosen[c1]-29*S.sqrt(lam)/(20160*S.pi**2))
    finite_ct = -4*S.pi*chosen[c0]*m/S.sqrt(lam)
    zero("volume_counterterm_finite_mass_piece", finite_ct-221*lam*m/(3780*S.pi))
    ren = S.simplify(outer(clock+Bct.subs(chosen)))
    zero("renormalized_local_clock_finite", ren
         -(lam*m*S.log(B*z)/(5*S.pi)+7*lam*m/(180*S.pi)))
    dm = S.symbols("dm")
    mass_bulk = outer(PfRN*(-2*dm/r)+QfRN*(2*dm/r**2))
    zero("normalizable_mass_bulk_boundary", mass_bulk+lam*dm/(36*S.pi))
    zero("normalizable_mass_W_boundary", outer(S.diff(BW,m))*dm+29*lam*dm/(945*S.pi))
    zero("normalizable_mass_CT_boundary", outer(S.diff(Bct.subs(chosen),m))*dm
         -221*lam*dm/(3780*S.pi))
    zero("normalizable_mass_boundary_cancels", mass_bulk
         +outer(S.diff(BW+Bct.subs(chosen),m))*dm)
    # Bulk value, not clock momentum: its vacuum is compensated by Lambda0.
    LW_RN = S.simplify(LW.subs(flat_alpha).subs(rn))
    bulk_series = outer(LW_RN)
    zero("W_bulk_cubic_action_divergence", bulk_series.coeff(z,-2)+13*lam**2/(240*S.pi))
    zero("W_bulk_linear_action_divergence", bulk_series.coeff(z,0)+lam/(60*S.pi))
    zero("W_action_cubic_divergence", -R(13,720)+R(463,7560)-R(653,15120))
    zero("vacuum_centered_total_cubic_divergence", R(653,15120)+R(11,720)-R(221,3780))
    # Einstein sign: derive its alpha boundary momentum before adding GHY/CT.
    R4 = -f2-3*a1*f1-2*f0*(a2+a1**2)-4*f1/r-4*f0*a1/r+2*(1-f0)/r**2
    Lam0 = S.symbols("Lam0")
    LEH = S.exp(a0)*r**2*(R4-2*Lam0)/(4*G)
    PEH = (S.diff(LEH,a1)-dr(S.diff(LEH,a2))).subs(flat_alpha)
    zero("Einstein_bulk_alpha_momentum", PEH+r**2*f1/(4*G))
    GHY = r**2*S.diff(F,r)/(4*G)+r*F/G
    CTEH = -S.sqrt(F)*(S.sqrt(lam)*r**2+1/(2*S.sqrt(lam)))/G
    EHclock = S.simplify(PEH.subs(rn)+GHY+CTEH)
    EEH = S.sqrt(F)*(r**2*S.sqrt(lam)+1/(2*S.sqrt(lam))-r*S.sqrt(F))/G
    zero("Einstein_clock_is_negative_energy", EHclock+EEH)
    zero("Einstein_positive_mass_sign", outer(EHclock)+m/G)
    dlambda = 11*G*lam**2/(360*S.pi)
    zero("vacuum_Einstein_bulk_excess", 3*dlambda/(2*G)-11*lam**2/(240*S.pi))


def rational_benchmark():
    k=S.Integer(2000)
    check("benchmark_X_lower_guard", (k**3/2)*(1+3*k**3/2)<3*k**6)
    check("benchmark_X_upper_guard", (2*k**3)*(1+6*k**3)>R(22,7)*k**6)
    check("vacuum_ratio_above_half", R(11,360*3)*k**-10<1)
    check("lambda_bclassical_square_above_one", k**3/4>1)
    check("logarithm_abs_bound", R(8,3)**24>2*k**3)
    check("charge_bracket_bound", 12+1+R(151,180)<14)
    check("fixed_G_relative_charge_bound", R(14,15)*k**-10<R(1,10**33))
    check("conditional_scope_recorded", "physical selection of a gravitational boundary ensemble" in CONTRACT["excluded"])


def hashes(data):
    return {name: hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


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
            mellin_and_distributions()
            local_variation()
            rational_benchmark()
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
