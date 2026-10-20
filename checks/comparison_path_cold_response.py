"""S156: closed-interval transfer and actual-path cold canonical response."""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "comparison-path-cold-response"
INPUTS = [
    "checks/comparison_path_cold_response.py",
    "notes/comparison-path-cold-response.md",
    "notes/comparison-path-interior-response.md",
    "evidence/comparison-path-interior-response.json",
    "notes/cold-response-from-real-endpoints.md",
    "evidence/cold-response-from-real-endpoints.json",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/divided-jet-dirac-response.md",
    "notes/interior-pressure-and-bolt-primitive.md",
]
PROTOCOL = {
    "digits": 70, "ellipse_rho": 256, "chebyshev_degrees": [2, 8, 16, 32],
    "scalar_parameters": ["0", ".3", "1"], "contour_nodes": [64, 128],
    "contour_radius": ".01", "contour_tolerance": "1e-35",
    "scalar_steps": ["1e-7", "5e-8"], "difference_tolerance": "1e-7",
    "refinement_ratio": ".35", "refinement_floor": "1e-40",
    "endpoint_differences": "second-order one-sided within [0,1]",
    "pressure_parameters": ["0", ".4", "1"], "pressure_cuts": ["exp(-4)", ".1", "1"],
    "pressure_gamma": ".03", "pressure_gap": "exp(-79)",
    "pressure_steps": ["1e-4", "5e-5"], "gauss_nodes": [32, 64],
    "pressure_normalization": "h L_h sqrt(c)", "quadrature_tolerance": "1e-25",
    "fixtures_are_quantum_sources": False,
}
ROWS, NUMERICS = [], []


def check(name, passed, detail=""):
    if name in {r["name"] for r in ROWS}:
        raise AssertionError("duplicate check: "+name)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessor_checks():
    for stem, count in [
        ("comparison-path-interior-response", 126),
        ("cold-response-from-real-endpoints", 152),
        ("fixed-chart-global-constraint-transfer", 99),
    ]:
        run = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))["runs"][-1]
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(r["passed"] for r in run["checks"])
              and len({r["name"] for r in run["checks"]}) == count)
        check(stem+"_current_inputs", run["unchanged_inputs"] and
              run["input_sha256_before"] == run["input_sha256_after"] and all(
                  sha256((ROOT/p).read_bytes()).hexdigest() == h
                  for p, h in run["input_sha256_after"].items()))


def exact_checks():
    rho = Q(256)
    lammax = Q(1, 2)+(rho+1/rho)/4
    check("ellipse_parameter_radius", lammax == Q("64.5009765625"))
    check("ellipse_coefficient_cone", lammax*Q(".00013") < Q(".008386") < Q(".011"))
    eta = Q("8.5e-33")
    K = 2*eta*(Q(2, 3)*10**17+9*10**11*24*320**3+10**17)
    check("complete_complex_difference_bound", 4*K < Q("5e-11"))
    check("real_endpoint_scalar_amplitude", 2*eta*10**23*Q("1e-6") == Q("1.7e-15"))
    check("real_endpoint_bolt_amplitude", 2*eta*10**24*Q("1e-6") == Q("1.7e-14"))
    t = s.symbols("t", nonnegative=True)
    for power in [2, 4]:
        numerator = s.expand(-s.Rational(1, 4)*(91+t/2)+s.Rational(power, 2))
        check(f"logarithmic_tail_monotone_{power}",
              s.Poly(numerator, t).all_coeffs() == [-s.Rational(1, 8), -s.Rational(91, 4)+s.Rational(power, 2)]
              and all(c < 0 for c in s.Poly(numerator, t).all_coeffs()))
    k = s.symbols("k", positive=True)
    for j in [1, 2]:
        radius = 2*j*j/k**2
        coefficient = s.simplify(2**j*s.factorial(j)/radius**j)
        target = k**2 if j == 1 else k**4/8
        check(f"lambda_cauchy_factor_{j}", s.simplify(coefficient-target) == 0)
    # The one-sided formulas cancel degree <= j+1 Taylor errors.
    e = s.symbols("e", positive=True)
    coeff = s.symbols("c0:5")
    polynomial = sum(coeff[i]*t**i for i in range(5))
    values = [polynomial.subs(t, n*e) for n in range(4)]
    d1 = (-3*values[0]+4*values[1]-values[2])/(2*e)
    d2 = (2*values[0]-5*values[1]+4*values[2]-values[3])/e**2
    check("one_sided_first_order_two", s.expand(d1-coeff[1]) == -2*coeff[3]*e**2-6*coeff[4]*e**3)
    check("one_sided_second_order_two", s.expand(d2-2*coeff[2]) == -22*coeff[4]*e**2)

    F1, F2 = Q(".000253"), Q("2.5e-11")
    inputs = [(Q("3e-7"), Q("3.5e-12"), Q("9e-10")),
              (Q("5e-4"), Q("3.5e-11"), Q("9e-9"))]
    bounds = [(Q("1.1e-6"), Q("1.7e-10"), Q("3.1e-9")),
              (Q(".0017"), Q("2.6e-7"), Q("3.1e-8"))]
    for tag, source, bound in zip(["cold", "bolt"], inputs, bounds):
        z0, z1, z2 = source
        actual = [Q("3.4")*z0, Q("3.4")*z1+2*F1*z0,
                  Q("3.4")*z2+4*F1*z1+2*F2*z0]
        for j in range(3):
            check(f"{tag}_profile_budget_{j}", actual[j] <= bound[j], str(actual[j]))
        check(f"{tag}_kernel_cross_terms_present", actual[2] > Q("3.4")*z2)
    check("quarter_power_bolt_integral",
          s.integrate(t**(-s.Rational(3, 4)), (t, 0, 1)) == 4)
    check("quarter_power_bolt_constant", 3*81*4 == 972)
    dg = Q("3e-7")
    N0, N1, N2 = bounds[0]
    B0, B1, B2 = bounds[1]
    m2 = Q("1.003")*N2+Q("4.02")*dg*N1+Q("5.54")*dg*dg*N0
    c2 = Q("1.006")*N2+Q("10.08")*dg*N1+Q("26.5")*dg*dg*N0
    b2 = Q("1.006")*B2+Q("10.08")*dg*B1+Q("26.5")*dg*dg*B0
    outputs = [46*m2/2, Q(".000002")*m2/2, (13122*c2+972*b2)/2]
    for tag, value, bound in zip(["mass", "density", "lapse"], outputs,
                                  [Q("7.2e-8"), Q("3.2e-15"), Q("3.6e-5")]):
        check("cold_constraint_"+tag, value < bound, str(value))
    check("geometry_cross_terms_retained", m2 > Q("1.003")*N2 and c2 > Q("1.006")*N2
          and b2 > Q("1.006")*B2)
    check("lapse_taylor_allowance_is_coarse", Q("3.6e-5")/2 > 1000*Q("1.2e-8"))
    NUMERICS.append({"cold_constraint_second_exact": [str(v) for v in outputs],
                     "full_goal_complete": False, "general_family_proved": False,
                     "homogeneous_density_included_in_deviation_bound": False})


def ellipse_checks():
    rho = mp.mpf(256)
    log_rho = mp.log(rho)
    check("ellipse_cauchy_disks_fit", log_rho > 4)
    check("cold_log_budget", 80+mp.log(mp.mpf("5e-11")/mp.mpf("1.7e-15")) < 91)
    check("bolt_log_budget", 80+mp.log(mp.mpf("5e-11")/mp.mpf("1.7e-14")) < 91)
    kk = 91/log_rho
    values = [mp.mpf("1.7e-15")*mp.exp(2)*kk**2,
              mp.mpf("1.7e-15")*mp.exp(4)*kk**4/8]
    for j, val, bound in zip([1, 2], values, ["3.5e-12", "9e-10"]):
        check(f"scalar_response_constant_{j}", val < mp.mpf(bound), mp.nstr(val, 35))
    # Independent foci-distance and inverse-Joukowski evaluations.
    for modulus in [mp.mpf("1.001"), mp.mpf(2), rho]:
        for index, theta in enumerate([mp.mpf(".13"), mp.mpf("1.2"), mp.mpf("2.8")]):
            w = modulus*mp.exp(1j*theta)
            z = (w+1/w)/2
            focus = (abs(z-1)+abs(z+1))/2
            check(f"ellipse_green_identity_{modulus}_{index}",
                  abs(mp.acosh(focus)-mp.log(modulus)) < mp.mpf("1e-60"))
    for center in [-1, 0, 1]:
        for j in [1, 2]:
            radius = 2*j*j/kk**2
            errors = []
            for n in range(32):
                z = center+radius*mp.exp(2j*mp.pi*(n+mp.mpf(".25"))/32)
                potential = mp.acosh((abs(z-1)+abs(z+1))/2)
                errors.append(potential-mp.sqrt(2*radius))
            check(f"closed_interval_cauchy_disk_{center}_{j}", max(errors) < 0)
    NUMERICS.append({"unrounded_scalar_response": [mp.nstr(v, 35) for v in values]})


def finite_difference(fun, x, h, derivative):
    if x == 0 or x == 1:
        step = h if x == 0 else -h
        assert all(0 <= x+n*step <= 1 for n in range(4))
        vals = [fun(x+n*step) for n in range(4)]
        if derivative == 1:
            return (-3*vals[0]+4*vals[1]-vals[2])/(2*step)
        return (2*vals[0]-5*vals[1]+4*vals[2]-vals[3])/step**2
    assert 0 <= x-h <= x+h <= 1
    if derivative == 1:
        return (fun(x+h)-fun(x-h))/(2*h)
    return (fun(x+h)-2*fun(x)+fun(x-h))/h**2


def difference_checks(tag, fun, x, target, derivative, steps):
    errors = [abs(finite_difference(fun, x, mp.mpf(h), derivative)-target) for h in steps]
    check(tag+"_difference", max(errors) < mp.mpf(PROTOCOL["difference_tolerance"])*(1+abs(target)),
          {"target": mp.nstr(target, 25), "errors": [mp.nstr(e, 12) for e in errors]})
    check(tag+"_refinement", errors[1] < mp.mpf(PROTOCOL["refinement_ratio"])*errors[0]
          or max(errors) < mp.mpf(PROTOCOL["refinement_floor"]))


def scalar_fixtures():
    x, w = s.symbols("x w")
    for degree in PROTOCOL["chebyshev_degrees"]:
        poly = s.chebyshevt(degree, 2*x-1).expand()
        identities = [s.diff(poly, x, j) for j in [1, 2]]
        check(f"chebyshev_endpoint_first_{degree}", identities[0].subs(x, 1) == 2*degree**2)
        check(f"chebyshev_endpoint_second_{degree}",
              identities[1].subs(x, 1) == s.Rational(4, 3)*degree**2*(degree**2-1))
        check(f"chebyshev_joukowski_{degree}", s.expand(
            s.chebyshevt(degree, (w+1/w)/2)-(w**degree+w**(-degree))/2) == 0)
        funcs = [s.lambdify(x, p, "mpmath") for p in identities]
        # Evaluation via recurrence is independent of polynomial coefficients.
        fun = lambda lam: mp.chebyt(degree, 2*lam-1)
        log_ratio = mp.log((mp.mpf(256)**degree+mp.mpf(256)**(-degree))/2)
        kk = max(1, log_ratio/mp.log(256))
        for parameter in PROTOCOL["scalar_parameters"]:
            xx = mp.mpf(parameter)
            for j in [1, 2]:
                target = funcs[j-1](xx)
                bound = mp.exp(2)*kk**2 if j == 1 else mp.exp(4)*kk**4/8
                tag = f"scalar_{degree}_{parameter}_{j}"
                check(tag+"_interval_bound", abs(target) <= bound)
                radius = mp.mpf(PROTOCOL["contour_radius"])
                for count in PROTOCOL["contour_nodes"]:
                    total = mp.mpc(0)
                    for n in range(count):
                        unit = mp.exp(2j*mp.pi*n/count)
                        total += fun(xx+radius*unit)/unit**j
                    contour = mp.factorial(j)*total/count/radius**j
                    check(tag+f"_contour_{count}", abs(contour-target) <
                          mp.mpf(PROTOCOL["contour_tolerance"])*(1+abs(target)))
                difference_checks(tag, fun, xx, target, j, PROTOCOL["scalar_steps"])
        NUMERICS.append({"degree": degree, "second_endpoint_over_n_squared":
                         str(Q(4, 3)*(degree**2-1)), "amplitude_normalized": True})
    check("endpoint_uniform_n_squared_bound_rejected",
          Q(4, 3)*(32**2-1) > 1000 and Q(4, 3)*(32**2-1) > 4*Q(4, 3)*(8**2-1))


def pressure_fixtures():
    nodes = {n: mp.gauss_quadrature(n, "legendre") for n in PROTOCOL["gauss_nodes"]}

    def gauss(fun, count):
        xs, ws = nodes[count]
        return sum(weight*fun((point+1)/2)/2 for point, weight in zip(xs, ws))

    gamma = mp.mpf(PROTOCOL["pressure_gamma"])
    gap = mp.exp(-79)
    Lh = 1+mp.log(mp.mpf(".01")/gap)
    cuts = [mp.exp(-4), mp.mpf(".1"), mp.mpf(1)]
    for cutlabel, cut in zip(PROTOCOL["pressure_cuts"], cuts):
        g = gamma*cut**2
        # Already normalized by gap*Lh*sqrt(cut); no tiny absolute tolerances.
        def pressure(lam, count):
            H = 1+mp.mpf(".2")*lam+mp.mpf(".1")*lam**2
            return H*gauss(lambda v: 2*v**4*mp.exp(2*lam*g*(v**4-1))
                          *(1+2*lam*g*v**4), count)
        for parameter in PROTOCOL["pressure_parameters"]:
            lam = mp.mpf(parameter)
            H = 1+mp.mpf(".2")*lam+mp.mpf(".1")*lam**2
            H1 = mp.mpf(".2")*(1+lam)
            H2 = mp.mpf(".2")

            def density(v, omit=False):
                d = g*(v**4-1)
                M = 1+2*lam*g*v**4
                bracket = M*H2
                if not omit:
                    bracket += (4*d*M+4*g*v**4)*H1+(4*d*d*M+8*d*g*v**4)*H
                return 2*v**4*mp.exp(2*lam*d)*bracket

            target = mp.quad(density, [0, mp.mpf(".5"), 1])
            tag = f"pressure_{cutlabel}_{parameter}"
            for count in PROTOCOL["gauss_nodes"]:
                check(tag+f"_second_quadrature_{count}",
                      abs(gauss(density, count)-target) < mp.mpf("1e-25"))
            check(tag+"_value_node_agreement",
                  abs(pressure(lam, 32)-pressure(lam, 64)) < mp.mpf("1e-25"))
            difference_checks(tag, lambda q: pressure(q, 64), lam, target, 2, PROTOCOL["pressure_steps"])
            omitted = mp.quad(lambda v: density(v, True), [0, mp.mpf(".5"), 1])
            check(tag+"_kernel_omission_detected", abs(target-omitted) > mp.mpf("1e-8"))
            NUMERICS.append({"pressure_parameter": parameter, "cut": cutlabel,
                             "physical_normalization": mp.nstr(gap*Lh*mp.sqrt(cut), 25),
                             "normalized_second_response": mp.nstr(target, 25)})
    c = s.symbols("c", positive=True)
    H = s.symbols("H", positive=True)
    raw = s.integrate(c**s.Rational(3, 2), c)/c**2*H
    null = s.simplify(2*raw-H*s.sqrt(c))
    check("zero_gamma_pressure_exact", s.simplify(raw-s.Rational(2, 5)*H*s.sqrt(c)) == 0)
    check("zero_gamma_null_exact", s.simplify(null+H*s.sqrt(c)/5) == 0)
    check("zero_gamma_absolute_bolt_integral",
          s.integrate(-null/c, (c, 0, 1)) == 2*H/5)


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = PROTOCOL["digits"]
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S156", "input_count": len(before)}), flush=True)
    try:
        predecessor_checks()
        exact_checks()
        ellipse_checks()
        scalar_fixtures()
        pressure_fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S156", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error, "scope": "actual-path cold canonical deviation response including both parameter endpoints",
        "unresolved": ["outer and other source allocations", "general family and full higher-order inverse",
                       "coupled solve", "physical completion requirements"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out = ROOT/"evidence"/(STEM+".json")
    history = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs": []}
    history["runs"].append(record)
    out.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
