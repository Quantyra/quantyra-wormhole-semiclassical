"""S159: exact LLL HH stress and reference/homogeneous allocations."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "comparison-path-lll-reference-response"
INPUTS = [
    "checks/comparison_path_lll_reference_response.py",
    "notes/comparison-path-lll-reference-response.md",
    "notes/comparison-path-interior-response.md",
    "evidence/comparison-path-interior-response.json",
    "notes/comparison-path-cold-response.md",
    "evidence/comparison-path-cold-response.json",
    "notes/comparison-path-outer-recombination.md",
    "evidence/comparison-path-outer-recombination.json",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-horizon-source-defect.md",
    "notes/stationary-reference-cold-bolt.md",
]
PROTOCOL = {
    "digits": 70, "chi": [".001", ".2", "1", "4"], "parameters": ["0", ".4", "1"],
    "gauss_nodes": [32, 64], "ward_tolerance": "1e-28",
    "ward_equal_panels": 8,
    "transformation_tolerance": "1e-35",
    "contour_nodes": [64, 128], "contour_radius": ".01",
    "real_steps": [".001", ".0005"], "second_tolerance": "1e-6",
    "refinement_ratio": ".35", "refinement_floor": "1e-35",
    "endpoint_differences": "second-order one-sided within [0,1]",
    "fixtures_are_physical_sources": False,
}
ROWS, NUMERICS = [], []
Q = s.Rational


def check(name, passed, detail=""):
    assert name not in {x["name"] for x in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessors():
    for stem, count in [("comparison-path-interior-response", 126),
                        ("comparison-path-cold-response", 242),
                        ("comparison-path-outer-recombination", 242),
                        ("fixed-chart-global-constraint-transfer", 99)]:
        run = json.loads((ROOT/"evidence"/(stem+".json")).read_text())["runs"][-1]
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(x["passed"] for x in run["checks"])
              and len({x["name"] for x in run["checks"]}) == count)
        check(stem+"_current_inputs", run["unchanged_inputs"]
              and run["input_sha256_before"] == run["input_sha256_after"]
              and all(sha256((ROOT/p).read_bytes()).hexdigest() == h
                      for p, h in run["input_sha256_after"].items()))


def exact_checks():
    a, a1, a2, beta, k = s.symbols("a a1 a2 beta k")
    pressure = beta*(1-a1*a1)/a**2
    curvature = -2*a2/a
    D = lambda f: s.diff(f, a)*a1+s.diff(f, a1)*a2
    check("HH_Ward_antiderivative", s.expand(D(a*a*pressure)-a*a1*beta*curvature) == 0)
    check("unique_regular_bolt_constant", (beta*(k-a1*a1)).subs(a1, 1) == beta*(k-1))
    w1, w2 = s.symbols("w1 w2")
    check("proper_density_identity",
          s.expand((pressure-beta*curvature).subs({a1: a*w1, a2: a*(w2+w1*w1)})
                   -beta*(a**-2+2*w2+w1*w1)) == 0)
    check("proper_null_identity",
          s.expand((2*pressure-beta*curvature).subs({a1: a*w1, a2: a*(w2+w1*w1)})
                   -2*beta*(a**-2+w2)) == 0)
    sig, sp, spp, coth, csch2 = s.symbols("sig sp spp coth csch2")
    wprime = s.exp(-sig)*(sp+coth)
    wsecond = s.exp(-2*sig)*(spp-csch2-sp*(sp+coth))
    ainv = s.exp(-2*sig)*csch2
    # coth^2-csch^2=1 is the only hyperbolic identity required.
    density = s.expand((ainv+2*wsecond+wprime*wprime)*s.exp(2*sig)).subs(csch2, coth*coth-1).expand()
    null = s.expand((ainv+wsecond)*s.exp(2*sig)).subs(csch2, coth*coth-1).expand()
    check("common_coordinate_density", density == 1+2*spp-sp*sp)
    check("common_coordinate_null", null == spp-coth*sp-sp*sp)
    u, u1, u2, lam = s.symbols("u u1 u2 lam")
    wpnew = s.exp(-lam*u)*(w1+lam*u1)
    wppnew = s.exp(-2*lam*u)*(w2+lam*u2-lam*u1*(w1+lam*u1))
    oldrho, oldnull = ainv+2*w2+w1*w1, 2*(ainv+w2)
    newrho = s.exp(-2*lam*u)*ainv+2*wppnew+wpnew**2
    newnull = 2*(s.exp(-2*lam*u)*ainv+wppnew)
    check("fixed_seed_density_transform", s.simplify(newrho-s.exp(-2*lam*u)*(oldrho+2*lam*u2-lam**2*u1**2)) == 0)
    check("fixed_seed_null_transform", s.simplify(newnull-s.exp(-2*lam*u)*(oldnull+2*lam*(u2-w1*u1)-2*lam**2*u1**2)) == 0)
    t, p0, pv, dp, K = s.symbols("t p0 pv dp K")
    f = -1+4/t-3/t**2
    pa = p0/t**4+pv*(1-1/t**4)+dp*(-1/t**4+4/t**5-3/t**6)
    ps = -p0/t**4+pv*(1+1/t**4)+dp*(1/t**4-6/t**5+6/t**6)
    ref = K*(p0-pv+pv*t**4+dp*f)
    check("reference_orbit_density", s.expand(K*t**4*pa-ref) == 0)
    check("reference_radius_Ward", s.simplify(ps-pa-t*s.diff(pa, t)/2) == 0)
    check("reference_product_orbit", pa.subs(t, 1) == p0)
    check("reference_product_radius", ps.subs(t, 1) == dp-p0+2*pv)
    check("nonvacuum_reference", s.expand(ref-K*pv*t**4-K*(p0-pv+dp*f)) == 0)
    rm, rl = s.symbols("rm rl")
    rstar = -K*p0-beta
    cold = rl-beta+K*(pv*(t**4-1)+dp*f)
    check("cold_homogeneous_cancellation", s.expand(cold-(rstar+rl+ref)) == 0)
    check("cold_total_equals_full", s.expand((rm-rstar)+cold-(rm+rl+ref)) == 0)
    check("middle_total_equals_full", s.expand(rm+rl+ref-(rm+rl+K*t**4*pa)) == 0)
    weyl = s.symbols("weyl")
    check("cold_middle_interface_with_Weyl",
          s.expand((rm-rstar)+cold+weyl-(rm+rl+ref+weyl)) == 0)
    check("outer_vacuum_not_counted_twice",
          s.expand((rm+weyl+K*pv*t**4)+(rl+K*(p0-pv+dp*f))
                   -(rm+weyl+rl+K*t**4*pa)) == 0)
    check("omitted_frozen_massive_density_detected", (beta+ref.subs(t, 1)) != 0)
    check("cold_center_zero", cold.subs({t: 1, rl: beta}) == 0)
    # Positive finite chart majorants, divided by the cold weight before evaluating.
    X, ax = s.symbols("X ax", positive=True)
    H = (1+Q("2e-6")*ax)*X+Q("3e-6")*ax*X**2+Q("2e-11")*ax**2*X**3
    unk = Q("1e-6")*ax*sum((100*H)**j for j in range(4))
    seed = Q(".026")*((1-4*H)**-2-(1-4*X)**-2)
    jac = sum((Q("6e-6")*ax*H)**j/j for j in range(1, 4))
    total = s.series(unk+seed+jac, X, 0, 4).removeO().expand()
    for j in [1, 2, 3]:
        normalized = s.cancel(total.coeff(X, j)*s.factorial(j)/ax)
        val = normalized.subs(ax, Q(".00331"))
        check(f"cold_fixed_chart_ordinary_jet_{j}", val < 7, str(val))
    check("cold_covariant_third_conversion",
          Q(".026")*3*4**2*2 < 3 and 7+120*Q("2e-6") < 10)


def budgets():
    eta, gamma = Q("8.5e-33"), Q("6.8e-28")
    beta = eta*400/36
    refscale = 8*Q(22, 7)*gamma
    check("massless_beta_guard", beta < Q("1e-31"))
    check("reference_scale_guard", refscale < Q("1.71e-26"))
    low = mp.mpf("1.2e-8")
    check("complex_frame_guard", mp.exp(66*low) < mp.mpf("1.000001"))
    check("complex_radius_guard", mp.exp(33*low)-1 < mp.mpf("4e-7"))
    check("cold_inverse_radius", 1/(1-mp.mpf("4e-7")) < mp.mpf("1.000001"))
    check("cold_radius_fourth_difference", 4*mp.mpf("1.002")**3*mp.mpf("1.000001") < mp.mpf("4.1"))
    check("cold_reference_shape_difference",
          (4*mp.mpf("1.000001")**2+6*mp.mpf("1.000001")**3)*mp.mpf("1.000001") < mp.mpf("10.1"))
    check("middle_radius_upper", 40*mp.exp(33*low) < 41)
    check("outer_inverse_radius", mp.exp(33*low)/40 < mp.mpf(".026"))
    check("outer_reference_shape", 1+4*Q(".026")+3*Q(".026")**2 < 2)
    check("bolt_seed_comparison", mp.exp(mp.mpf(".05"))*mp.sqrt(mp.mpf("1.01004")) < mp.mpf("1.1"))
    check("bolt_quadratic_seed_modulus", mp.cosh(mp.mpf("1.1"))-1 < 2)
    check("seed_null_bolt_modulus", 34*2 < 100)
    check("cold_curvature_seed_deviation", 17+6*Q("1.001")*Q("1e-30") < 18)
    yy = mp.mpf(".081")
    LL = 1-mp.log(yy)
    check("outer_direction_radial_Hessian",
          6*mp.mpf("1e-11")*yy+mp.mpf("2e-8")*2*100**2*yy**2*LL < mp.mpf("1e-4"))
    check("outer_direction_tangential_Hessian",
          yy*(3*mp.mpf("1e-11")*yy**2+mp.mpf("2e-8")*100*yy**3*LL) < mp.mpf("1e-4"))
    check("middle_fixed_normal_Hessian", Q("1.3e-4")*2*100**2 < 10)
    ax = Q(".00331")
    # Formula (3): product rescaling, changed Hessian and gradient-square.
    rho_bound = Q("1.000001")*(27+66*Q("1e-6")+2*33*10+33**2*100*ax)
    null_bound = Q("1.000001")*(36+2*33*20+2*33**2*100*ax)
    check("cold_massless_density_complex", rho_bound < 10**5, str(rho_bound))
    check("cold_massless_null_complex", null_bound < 10**5, str(null_bound))
    bolt = Q("1.000001")*(100+2*33*44+2*33**2*200*2*ax*Q("1.1"))
    check("bolt_massless_complex", bolt < 10**5, str(bolt))
    middle = Q("1.000001")*(4+2*33*20+2*33**2*Q("1e-4")**2)
    outer = Q("1.000001")*(20+2*33*Q("2e-4")+2*33**2*Q("1e-4")**2)
    check("middle_massless_complex", beta*middle < Q("1e-26"))
    check("outer_massless_complex", beta*outer < Q("1e-26"))
    check("cold_full_allocation", 10**5*beta+refscale*(Q("4.1")+3*Q("10.1")) < Q("1e-23"))
    check("bolt_full_allocation", 10**5*beta < Q("1e-24"))
    check("middle_full_allocation",
          Q("1e-26")+refscale*(3301+41**4+3*(1+4*Q("1.001")+3*Q("1.001")**2)) < Q("1e-18"))
    check("outer_full_allocation", Q("1e-26")+refscale*(3301+6) < Q("1e-21"))
    n0 = list(map(Q, ["1e-23", "1e-24", "1e-18", "1e-21"]))
    n1 = [v*Q(".032") for v in n0]
    n2 = [v*Q(".002") for v in n0]
    for i in range(4):
        check(f"first_Cauchy_row_{i}", n0[i]/32 < n1[i])
        check(f"second_Cauchy_row_{i}", 2*n0[i]/32**2 < n2[i])
    dg = Q("3e-7")
    mass1 = [Q("1.003")*b+Q("2.01")*dg*a for a, b in zip(n0, n1)]
    lapse1 = [Q("1.006")*b+Q("5.04")*dg*a for a, b in zip(n0, n1)]
    mass2 = [Q("1.003")*c+Q("4.02")*dg*b+Q("5.54")*dg**2*a for a, b, c in zip(n0, n1, n2)]
    lapse2 = [Q("1.006")*c+Q("10.08")*dg*b+Q("26.5")*dg**2*a for a, b, c in zip(n0, n1, n2)]
    matrix = [[46, 0, 279, Q(".00016")], [Q(".000002"), 0, 1, Q(".835")],
              [13122, 486, 1300, Q(".00003")]]
    for j, mass, lapse, limits in [(1, mass1, lapse1, ["4.6e-18", "1.7e-20", "2.2e-17"]),
                                   (2, mass2, lapse2, ["2.9e-19", "1.1e-21", "1.4e-18"])]:
        vals = []
        for i, tag in enumerate(["mass", "density", "lapse"]):
            value = sum(a*b for a, b in zip(matrix[i], lapse if i == 2 else mass))/2
            check(f"constraint_{tag}_response_{j}", value < Q(limits[i]), str(value))
            vals.append(str(value))
        NUMERICS.append({"constraint_order": j, "exact_bounds": vals})


def fixtures():
    chi, lam = s.symbols("chi lam", positive=True)
    sig0 = Q(".02")*chi**2/(1+chi**2)
    u = Q(".003")*chi**2/(1+chi**2)+Q(".002")/(1+chi**2)
    sig = sig0+lam*u
    a = s.exp(sig)*s.sinh(chi)
    D = lambda f: s.exp(-sig)*s.diff(f, chi)
    ap, app = D(a), D(D(a))
    trace = -2*app/a
    P = (1-ap**2)/a**2
    direct = [P-trace, 2*P-trace]
    rho0 = s.exp(-2*sig0)*(1+2*s.diff(sig0, chi, 2)-s.diff(sig0, chi)**2)
    null0 = 2*s.exp(-2*sig0)*(s.diff(sig0, chi, 2)-s.coth(chi)*s.diff(sig0, chi)-s.diff(sig0, chi)**2)
    Du = s.exp(-sig0)*s.diff(u, chi)
    D2u = s.exp(-2*sig0)*(s.diff(u, chi, 2)-s.diff(sig0, chi)*s.diff(u, chi))
    Dw0 = s.exp(-sig0)*(s.diff(sig0, chi)+s.coth(chi))
    transformed = [s.exp(-2*lam*u)*(rho0+2*lam*D2u-lam**2*Du**2),
                   s.exp(-2*lam*u)*(null0+2*lam*(D2u-Dw0*Du)-2*lam**2*Du**2)]
    # Direct common-coordinate formula, independent of transformed seed data.
    rho = s.exp(-2*sig)*(1+2*s.diff(sig, chi, 2)-s.diff(sig, chi)**2)
    null = 2*s.exp(-2*sig)*(s.diff(sig, chi, 2)-s.coth(chi)*s.diff(sig, chi)-s.diff(sig, chi)**2)
    bolt_series = s.series(null, chi, 0, 4).removeO().expand()
    check("fixture_bolt_null_constant_zero", bolt_series.coeff(chi, 0) == 0)
    check("fixture_bolt_null_linear_zero", bolt_series.coeff(chi, 1) == 0)
    check("fixture_bolt_quadratic_nonzero", bolt_series.coeff(chi, 2) != 0)
    free_state = 1/a**2
    check("free_state_pole_detected", s.simplify(s.limit(chi**2*free_state, chi, 0)-s.exp(-Q(".004")*lam)) == 0)
    t0 = 1+Q(".01")*s.cosh(chi)
    V = Q(".004")/(1+chi**2)
    t = 1+(t0-1)*s.exp(lam*V)
    pv, dp, p0 = -Q(".2"), Q(".03"), s.Integer(7)
    f = -1+4/t-3/t**2
    ref = p0-pv+pv*t**4+dp*f
    profiles = [rho-1+pv*(t**4-1)+dp*f, rho+ref,
                rho+p0-pv+dp*f, null]
    names = ["cold_density", "middle_density", "outer_density", "null"]
    value_fun = s.lambdify((chi, lam), profiles, "mpmath", cse=True)
    second_fun = s.lambdify((chi, lam), [s.diff(f, lam, 2) for f in profiles], "mpmath", cse=True)
    direct_fun = s.lambdify((chi, lam), direct, "mpmath", cse=True)
    transform_fun = s.lambdify((chi, lam), transformed, "mpmath", cse=True)
    pressure_fun = s.lambdify((chi, lam), P, "mpmath", cse=True)
    a2_fun = s.lambdify((chi, lam), a*a, "mpmath", cse=True)
    kernel = s.diff(a*a, chi)*trace
    kernel_fun = s.lambdify((chi, lam), kernel, "mpmath", cse=True)
    nodes = {n: mp.gauss_quadrature(n, "legendre") for n in PROTOCOL["gauss_nodes"]}
    def finite_second(fun, q, h):
        if q == 0:
            return (2*fun(q)-5*fun(q+h)+4*fun(q+2*h)-fun(q+3*h))/h**2
        if q == 1:
            return (2*fun(q)-5*fun(q-h)+4*fun(q-2*h)-fun(q-3*h))/h**2
        return (fun(q+h)-2*fun(q)+fun(q-h))/h**2
    for clabel in PROTOCOL["chi"]:
        c = mp.mpf(clabel)
        for qlabel in PROTOCOL["parameters"]:
            q = mp.mpf(qlabel)
            tag = f"fixture_{clabel}_{qlabel}"
            p = pressure_fun(c, q)
            pad = mp.quad(lambda v: kernel_fun(v, q) if v else 0, [0, c/2, c])/(2*a2_fun(c, q))
            check(tag+"_Ward_adaptive", abs(p-pad) < mp.mpf("1e-28"), mp.nstr(abs(p-pad), 25))
            for n in PROTOCOL["gauss_nodes"]:
                xs, ws = nodes[n]
                panels = PROTOCOL["ward_equal_panels"]
                pg = sum(w*c/(2*panels)*kernel_fun(c*(panel+(z+1)/2)/panels, q)
                         for panel in range(panels) for z, w in zip(xs, ws))/(2*a2_fun(c, q))
                check(tag+f"_Ward_gauss_{n}", abs(p-pg) < mp.mpf("1e-28"), mp.nstr(abs(p-pg), 25))
            dval, tval = direct_fun(c, q), transform_fun(c, q)
            for i, label in enumerate(["density", "null"]):
                check(tag+"_"+label+"_fixed_seed_transform", abs(dval[i]-tval[i]) < mp.mpf("1e-35"))
            target = second_fun(c, q)
            radius = mp.mpf(".01")
            for n in PROTOCOL["contour_nodes"]:
                sums = [mp.mpc(0) for _ in profiles]
                for k in range(n):
                    phase = mp.exp(2j*mp.pi*k/n)
                    vals = value_fun(c, q+radius*phase)
                    for i, val in enumerate(vals):
                        sums[i] += val/phase**2
                for i, label in enumerate(names):
                    error = abs(2*sums[i]/n/radius**2-target[i])
                    check(tag+f"_{label}_contour_{n}", error < mp.mpf("1e-6")*(1+abs(target[i])),
                          mp.nstr(error, 25))
            for i, label in enumerate(names):
                errors = []
                for step in PROTOCOL["real_steps"]:
                    val = finite_second(lambda z: value_fun(c, z)[i], q, mp.mpf(step))
                    error = abs(val-target[i])
                    errors.append(error)
                    check(tag+f"_{label}_difference_{step}", error < mp.mpf("1e-6")*(1+abs(target[i])),
                          mp.nstr(error, 25))
                check(tag+f"_{label}_refinement", errors[1] < mp.mpf(".35")*errors[0]
                      or max(errors) < mp.mpf("1e-35"))
            NUMERICS.append({"chi": clabel, "lambda": qlabel,
                             "second_responses": [mp.nstr(x, 35) for x in target]})


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = 70
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S159", "input_count": len(before)}), flush=True)
    try:
        predecessors()
        exact_checks()
        budgets()
        fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S159", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error, "scope": "actual-path LLL/reference/frozen homogeneous density and null response",
        "unresolved": ["remaining Weyl and full radius allocations", "E0 and full assembly",
                       "general family and full inverse", "coupled solution and physical completion"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out = ROOT/"evidence"/(STEM+".json")
    history = json.loads(out.read_text()) if out.exists() else {"runs": []}
    history["runs"].append(record)
    out.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
