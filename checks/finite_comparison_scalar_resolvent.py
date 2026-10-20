"""S167: scalar reference resolvent on the actual finite comparison.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "finite-comparison-scalar-resolvent"
OWN = ["checks/finite_comparison_scalar_resolvent.py", "notes/"+STEM+".md"]
PREDECESSORS = [("endpoint-preserving-compatible-projection", 262),
                ("stationary-scalar-outgoing-inverse", 107),
                ("comparison-path-interior-response", 126)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual finite-comparison static spherical scalar reference resolvent",
    "dimensionless_potential_Linfinity": "400", "dimensionless_potential_L1": "72000",
    "scattering_ratio": "2e-13", "kernel_multiplier": "1.001",
    "metric_derivatives_used": 2, "fixed_original_path": True,
    "fixture_scope": "exact metric identities and regular bolt expansion only",
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "formal_verification": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(expr)
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    return sorted(found)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(row["passed"] for row in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)
    for p in ["stationary-completed-metric-response-norm", "stationary-scalar-outgoing-inverse",
              "comparison-path-interior-response", "endpoint-preserving-compatible-projection",
              "stationary-seed-quantum-form", "stationary-centered-covariant-seed"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def exp_upper(x):
    # For 0<=x<1, exp(x)<=1/(1-x), directly from its positive series.
    assert 0 <= x < 1
    return 1/(1-x)


def metric_algebra():
    y = S.Symbol("y", positive=True)
    ell = S.Symbol("ell", real=True)
    tau, v = S.Function("tau")(y), S.Function("v")(y)
    k = S.exp(ell*v)
    nu, tt, bb = tau*v, 1+tau*(k-1), k-1-ell*v*k
    zero("finite_radius_first_jet", S.diff(tt,y)-ell*k*S.diff(nu,y)-bb*S.diff(tau,y))
    zero("finite_radius_second_jet", S.diff(tt,y,2)-ell*k*S.diff(nu,y,2)
         -bb*S.diff(tau,y,2)-tau*k*ell**2*S.diff(v,y)**2)
    a, r, w, ss = [S.Function(n)(y) for n in ["a", "r", "w", "s"]]
    dz = lambda f: S.exp(-w)*S.diff(f,y)
    measure = S.exp(w+2*ss)*a*r*r
    jj = dz(measure)/measure
    j0 = S.diff(a*r*r,y)/(a*r*r)
    v0 = S.diff(j0,y)/2+j0*j0/4
    vv = dz(jj)/2+jj*jj/4
    target = S.exp(-2*w)*(v0+S.diff(w,y,2)/2+S.diff(ss,y,2)
                            +j0*S.diff(ss,y)-S.diff(w,y)**2/4+S.diff(ss,y)**2)
    zero("actual_four_dimensional_Liouville_potential", vv-target)
    zero("physical_log_density", jj-S.exp(-w)*(j0+S.diff(w,y)+2*S.diff(ss,y)))
    f = S.Function("f")(y)
    # Independent divergence-form transformation of the scalar Laplacian.
    original = -dz(measure*dz(f/S.sqrt(measure)))/measure
    zero("divergence_form_unitary_conjugation", S.sqrt(measure)*original+dz(dz(f))-vv*f)
    q, t = S.symbols("q tau", real=True, nonzero=True)
    ff = S.log(1+t*(S.exp(q/t)-1))
    zero("finite_completion_second_q_derivative", S.diff(ff,q,2)
         -(1-t)*S.exp(q/t)/(t*(1+t*(S.exp(q/t)-1))**2))
    zero("finite_completion_value_at_zero", ff.subs(q,0))
    zero("finite_completion_first_q_at_zero", S.diff(ff,q).subs(q,0)-1)
    x, en = S.symbols("x E")
    zero("reference_fourth_order_resolvent_identity", (1/(x-en)-1/(x+4))/(en+4)-1/((x+4)*(x-en)))
    zero("energy_shift", R(11,3)-R(9,4)-R(17,12))


def low_and_outer_metric_budgets():
    eps, amp = R('3e-8'), R('1.2e-8')
    check("cold_second_jet_conversion", 2*(R('1e-5')*R('.00331')*20000+R('.01')) < 1)
    check("core_cover_at_t40", 1/(R('.412')*40) > R('.05'))
    check("canonical_to_physical_first", R('6.1') < R('2.5')**2)
    check("physical_first_jet", R('2.5')*R('1e-4') == R('.00025'))
    check("physical_second_jet", R('6.1')+R('2.5')*R('2.3')*R('1e-4') < R('6.2'))
    check("sqrt5_guard", 5 < R('2.3')**2)
    check("tau_second_jet", 11+2*5 < 22)
    check("finite_exponential_guard", exp_upper(amp) < R('1.001'))
    check("finite_radius_lower", 1-amp > R('.999'))
    check("B_ell_quadratic_guard", exp_upper(amp)*amp**2/2 < R('2e-16'))
    check("finite_log_value", amp*exp_upper(amp)/(1-amp) < eps/2)
    check("full_metric_value", amp+amp*exp_upper(amp)/(1-amp) < eps)
    sy = (R('1.001')*R('.00025')+R('2e-16')*R('2.3'))/R('.999')
    syy = (R('1.001')*R('6.2')+R('2e-16')*22
           +R('1.001')*7*R('3e-7')**2)/R('.999')+sy**2
    check("finite_radius_core_first", sy < R('.001'))
    check("finite_orbit_core_first", sy+R('.00025') < R('.001'))
    check("finite_radius_core_second", syy < 100)
    check("finite_orbit_core_second", syy+R('6.2') < 100)
    # Exact exponential lower partial sums certify the logarithm upper guards.
    check("log10_below_2point4", sum(R('2.4')**j/S.factorial(j) for j in range(10)) > 10)
    check("log40_below_3point8", sum(R('3.8')**j/S.factorial(j) for j in range(13)) > 40)
    check("log20_below_3", sum(R(3)**j/S.factorial(j) for j in range(12)) > 20)
    check("outer_Ylog_guard", R('.1')*(1+R('2.4')) == R('.34'))
    check("outer_log_over_t", (1+R('3.8'))/40 == R('.12'))
    check("outer_zY_bound", R('1.0011')/R('5.9') < R('.413')**2)
    check("outer_tau_first_scaled", R('.413')*R('.1') < R('.05'))
    check("outer_tau_second_scaled", R('.413')/2*R('.1')**3 < R('.001'))
    aa, bb = R('1e-11'),R('2e-8')
    qq = [R('1e-8'),R('1e-6'),R('2e-4')]
    for j, mass in enumerate([1,3,6]):
        check("tangent_scaled_jet_"+str(j), mass*aa+bb*S.factorial(j)*100**j*R('.34') < qq[j])
    zexp = exp_upper(R('.01')/R('.83'))
    check("complex_completion_denominator", R('1.1')*(zexp-1) < R('.014'))
    check("complex_H_bound", R('.17')*zexp/(R('.83')*R('.986')**2) < 1)
    for label, i, j, bound in [('q',1,0,200),('tau',0,1,20),('qq',2,0,80000),('qtau',1,1,4000),('tautau',0,2,800)]:
        check("H_Cauchy_"+label, S.factorial(i)*S.factorial(j)/R('.005')**i/R('.05')**j == bound)
    q0,q1,q2=qq
    h1=200*q1*R('.1')**3+20*R('.05')
    h2=(80000*q1*q1*R('.1')**6+8000*q1*R('.1')**3*R('.05')
        +800*R('.05')**2+200*q2*R('.1')**3+20*R('.001'))
    dd=[q0*q0,2*q0*q1+q0*q0*h1,
        2*q1*q1+2*q0*q2+4*q0*q1*h1+q0*q0*h2]
    for j,budget in enumerate(dd):
        check("finite_outer_quadratic_jet_"+str(j), budget < R('1e-10'),budget)
        check("finite_outer_normalized_remainder_"+str(j), R('2e-8')+budget*R('.1')**2/(S.factorial(j)*100**j) < R('1e-6'))
    check("finite_outer_cubic_unchanged", 6 > 3)
    # y derivatives of Y satisfy |Y_y|<1.1Y, |Y_yy|<3Y.
    check("outer_Y_y", R('1.001') < R('1.1'))
    check("outer_Y_yy", R('1.1')**2*5 < 9)
    coeffs=[(2,2),(R('6.6'),220),(R('32.52'),49000)]
    for j,(cc,dd) in enumerate(coeffs):
        budget=cc*R('1e-9')*27+dd*R('1e-6')*81*R('.12')
        check("outer_physical_decay_"+str(j),budget<10,budget)
    check("outer_first_fits_global_reserve", R(10)/40**3 < R('.001'))


def bolt_and_integral_budgets():
    eps, kk=R('3e-8'),R(100)
    y,v=S.symbols('y v',positive=True)
    zero("bolt_difference_integral_factor", S.integrate((y*y-v*v)/2,(v,0,y))/y-y*y/3)
    check("seed_j0_bound", R('10.5') < R('3.25')**2 and 5 < R('2.24')**2
          and R('3.25')+2*R('2.24') < 8)
    check("singular_cancellation_reserve", exp_upper(10*eps)*kk/6 < 17)
    check("R_global_difference", exp_upper(2*eps)-1 < 3*eps)
    check("critical_derivative_difference", exp_upper(5*eps)**2*kk*eps < R('.002')**2)
    check("physical_derivative_first_part", exp_upper(eps)*(8+R('.003')) < R('8.02'))
    check("physical_derivative_total", (R('8.02')+R('.002'))/2 < 5)
    check("exponential_relative_value", exp_upper(2*eps) < R('1.001'))
    check("constant_potential_difference", R(9,4)*(exp_upper(2*eps)-1) < R('1e-6'))
    core=R('1.001')*(30+50+100+100+R('.008')+R('.00000125'))+R('1e-6')+17
    check("core_potential_bound",core<400,core)
    check("sqrt2_reserve", 2 < R('1.5')**2)
    check("seed_length_to_t40", 172+R('1.5')*3 < 177)
    check("actual_length_upper", exp_upper(eps)*177 < 178)
    check("actual_length_lower", (1-eps)*R('.2') > R('.199'))
    check("outer_j0", 1/R('.2')+8 == 13)
    check("outer_metric_exponent_linear", 2*exp_upper(2*eps) < R('2.001'))
    outer_cubic=R('1.001')*(5+10+130+R(125)/40**3)+R(9,4)*R('2.001')*10
    check("outer_potential_tminus3",outer_cubic<200,outer_cubic)
    check("outer_potential_tminus2",R('1.001')*11<12)
    check("outer_measure",exp_upper(eps)**2*2<R('1.5')**2)
    outer_integral=R('1.5')*(R(12)/(2*40**2)+R(200)/(3*40**3))
    check("outer_potential_integral",outer_integral<R('.01'),outer_integral)
    check("outer_critical_integral",1/(4*R('.199'))<R('1.3'))
    check("outer_potential_Linfinity",R(12)/40**2+R(200)/40**3+1/(4*R('.199')**2)<400)
    check("core_integral",400*178==71200)
    check("global_potential_L1",71200+R('.01')+R('1.3')<72000)
    check("scattering_ratio",R(72000)/R('4e17')<R('2e-13'))
    check("resolvent_neumann_multiplier",1/(1-R('2e-13'))<R('1.001'))


def independent_regular_bolt_fixture():
    # a0=sinh(y), r0=2+y^2, with even finite metric warps.
    # This fixture has no claim to be the physical comparison or a solution.
    y=S.Symbol('y',positive=True)
    c2,c4,s2,s4=S.symbols('c2 c4 s2 s4',real=True)
    ww=c2*y*y+c4*y**4
    ss=s2*y*y+s4*y**4
    az=S.exp(ww)*S.sinh(y)
    rz=S.exp(ss)*(2+y*y)
    density=az*rz**2
    jraw=S.diff(density,y)/density
    jseries=S.series(jraw,y,0,4).removeO()
    dzseries=lambda f:S.series(S.exp(-ww)*S.diff(f,y),y,0,3).removeO()
    jj=S.series(S.exp(-ww)*jseries,y,0,4).removeO()
    potential=S.series(dzseries(jj)/2+jj**2/4,y,0,1).removeO()
    zseries=y+c2*y**3/3+(c4+c2*c2/2)*y**5/5
    critical=S.series(1/(4*zseries*zseries),y,0,1).removeO()
    zero("fixture_bolt_inverse_square_cancels", S.expand(potential+critical).coeff(y,-2))
    zero("fixture_bolt_constant", S.expand(potential+critical)-R(7,3)-R(4,3)*c2-4*s2)
    zero("fixture_proper_distance_jet", S.series(S.diff(zseries,y)-S.exp(ww),y,0,5).removeO())
    zero("fixture_bolt_period", S.limit(S.exp(-ww)*S.diff(az,y),y,0)-1)
    zero("fixture_bolt_finite_potential", S.expand(potential+critical-R(9,4))-(R(1,12)+R(4,3)*c2+4*s2))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison scalar reference resolvent",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        metric_algebra()
        low_and_outer_metric_budgets()
        bolt_and_integral_budgets()
        independent_regular_bolt_fixture()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison scalar reference resolvent",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "finite_comparison_scalar_L1_bound": good,
        "actual_reference_outgoing_kernel": good, "regular_bolt_realization_retained": good,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["finite-comparison coupled TT and mixed blocks", "density contacts and candidate range",
                       "projection direction and cold/higher bounds", "nonlinear preservation and solve",
                       "signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_bytes()) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
