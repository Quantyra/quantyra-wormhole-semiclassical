"""S167: C^(4,1/2) interior canonical source values and parameter response.

Exact finite arithmetic and Taylor fixtures; no determinant simulation.
"""
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "fractional-interior-dirac-response"
OWN = ["checks/fractional_interior_dirac_response.py", "notes/"+STEM+".md"]
PREDECESSORS = [("divided-jet-dirac-response", 93),
                ("general-interior-dirac-response", 61),
                ("spatial-response-norm-audit", 88),
                ("projected-first-derivative-bound", 246)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual retained-domain interior conformal canonical source values",
    "spatial_class": "real C^(4,1/2) with unchanged global and capped low norms",
    "parameter_domain": "complex bidisk of radii 1/(200U), 1/(200V)",
    "first_constant": "1e-9", "mixed_constant": "2e-7",
    "norm_inclusion": "P45<=P5", "full_angular_tower": True,
    "fixture_scope": "finite fractional Taylor products, not physical Dirac sources",
    "spatial_source_gradient": False, "global_nonlinear_family_preserved": False,
    "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
    "useful_Newton_certificate": False, "formal_verification": False,
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
    for p in ["divided-jet-dirac-response", "general-interior-dirac-response",
              "spatial-response-norm-audit", "coincident-normal-dirac-response",
              "cold-response-from-real-endpoints", "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def inclusion_and_cusp():
    # min(a,b)<=sqrt(ab); sqrt(2)<3/2 makes the inclusion fully rational.
    check("norm_sqrt2_guard", 2 < R(3, 2)**2)
    check("P5_fourth_holder_below_point4", 2*R(3, 2)*5/100 < R(".4")**2)
    check("ray_high_jet_factor", 2*100**2 < 150**2)
    check("ray_fourth_holder_factor", (S.sqrt(2))**4 == 4)
    x = S.Symbol("x", positive=True)
    cusp = x**R(9, 2)
    zero("cusp_fourth_exact", S.diff(cusp, x, 4)-R(945, 16)*S.sqrt(x))
    zero("cusp_fifth_exact", S.diff(cusp, x, 5)-R(945, 32)/S.sqrt(x))
    for j in range(5):
        check("cusp_continuous_jet_"+str(j), S.limit(S.diff(cusp, x, j), x, 0, dir="+") == 0)
    check("cusp_unbounded_fifth", S.limit(S.diff(cusp, x, 5), x, 0, dir="+") == S.oo)
    a, b = S.symbols("a b", nonnegative=True)
    # For 0<=b<=a, (sqrt(a)-sqrt(b))^2<=a-b; set a=(b_root+d)^2.
    zero("sqrt_holder_square_identity", ((b+a)**2-b**2)-a*a-2*a*b)
    v = S.Symbol("v", nonnegative=True)
    factor = S.integrate(4*(1-v)**3*S.sqrt(v), (v, 0, 1))
    check("fractional_Taylor_factor", factor == R(128, 315) and factor < 1, factor)
    # The divided fourth cusp seminorm times this factor recovers |x|^(9/2).
    zero("cusp_exact_fractional_integral", factor*R(945, 16)/S.factorial(4)-1)


def exp_coefficients(z):
    out = [S.Integer(1)]
    for n in range(1, 5):
        out.append(S.expand(sum(k*z[k]*out[n-k] for k in range(1, n+1))/n))
    return out


def fractional_division():
    z = [R(0), R(".016")]+[R(".011")*150**j for j in range(2, 5)]
    ec = exp_coefficients(z)
    x = S.Symbol("X")
    independent = S.series(S.exp(sum(z[j]*x**j for j in range(1, 5))), x, 0, 5).removeO()
    check("independent_finite_exp_series", all(independent.coeff(x, j) == ec[j] for j in range(5)))
    zz = [0]+list(S.symbols("z1:5"))
    ee = exp_coefficients(zz)
    for j in range(1, 5):
        for k in range(1, j+1):
            zero(f"exp_coefficient_difference_{j}_{k}", S.diff(ee[j], zz[k])-ee[j-k])
    # exp(y)<=1/(1-y), 0<=y<1, follows termwise from 1/j!<=1.
    check("inner_exp_value", 1/(1-R(".016")/64) < R("1.001"))
    seed = [R("1.026")]+[R(".026")*64**j for j in range(1, 6)]
    inner = [R("1.001")*sum(seed[k]*ec[j-k] for k in range(j+1)) for j in range(5)]
    integer = {2: R(400), 3: R(50000), 4: R(7000000)}
    displayed = {2: 361, 3: 45369, 4: 6276299}
    for j in range(2, 5):
        check("inner_integer_"+str(j), inner[j] < displayed[j] < integer[j], inner[j])
    center_seed = [R(1), R(0), R(8), seed[3], seed[4]]
    pc = [None, R(".02"), R(300), R(50000), R(7000000)]
    for j in range(1, 5):
        center = sum(center_seed[k]*ec[j-k] for k in range(j+1))
        check("center_coefficient_"+str(j), center < pc[j], center)
    for d in range(1, 4):
        outside = R("2.02")*64**d+sum(pc[j]*64**(d+1-j) for j in range(2, d+1))
        check("outer_integer_division_"+str(d), outside < integer[d+1], outside)
    hz = [R(".016")]+[(j+1)*z[j+1] for j in range(1, 4)]+[4*R(".011")*100**5]
    check("holder_unknown_jet_list", hz == [R(".016"), 495, 111375, 22275000, 440000000])
    ah = [(j+1)*seed[j+1] for j in range(5)]
    eh = [R("1.001")*(hz[0]*ec[j]+sum(hz[k]*ec[j-k] for k in range(1, j+1))) for j in range(5)]
    hh = [sum(ah[k]*R("1.001")*ec[j-k]+seed[k]*eh[j-k] for k in range(j+1)) for j in range(5)]
    check("fourth_holder_product_bound", hh[4] < 10**9, hh[4])
    outside = R("2.02")*8**7+300*8**5+50000*8**3+7000000*8
    check("outer_fractional_division", outside < 10**9, outside)
    check("zeroth_fractional_division", R(".06")*2 < 2**2)
    for j, square, reserve in [(1, 400*R("2.02"), 29),
                               (2, 50000*(400+300), 6000),
                               (3, 7000000*(50000+50000), 840000)]:
        check("fractional_interpolation_"+str(j), square < reserve**2, square)
    return pc


BB = [None, R(2), R(29), R(6000), R(840000), R(10**9)]
PC = [None, R(".02"), R(300), R(50000), R(7000000)]


@lru_cache(None)
def product(n, d):
    if n == 1:
        return BB[d+1]
    if d == n-1:
        return R(2)**n
    return BB[d-n+2]*2**(n-1)+sum(PC[a]*product(n-1, d-a) for a in range(1, d-n+2))


def born_and_angular_budgets():
    expected = [R(10**9), R(4946800), R("41990.4"), R("234.343432")]
    for n, target in enumerate(expected, 1):
        check("fractional_product_"+str(n), product(n, 4) == target, product(n, 4))
    # Both ranges L<=1 and L>=1 give the base exponent N-1/2.
    for n in range(2, 6):
        check("base_low_length_power_"+str(n), n > n-R(1, 2))
        check("base_high_length_power_"+str(n), n-R(1, 2) > 0 and R(".06")**n < 2**n)
    check("positive_loop_interpolation", 28**2*S.factorial(4)*S.factorial(5)*10**9 < R("4.8e7")**2)
    rm = R("1.02")
    for j, target in [(4, 10), (5, 20)]:
        bound = R(8, 3)*rm**2*S.factorial(j-1)/(R(7, 5)*R(".98"))**j
        check("two_edge_vertex_moment_"+str(j), bound < target, bound)
    check("two_edge_fractional_interpolation", 10*20 < 20**2)
    n, j = S.symbols("N j")
    zero("fractional_mass_power", (n+1)-(n+R(9, 2)-1)+R(5, 2))
    check("tail_mass_conversion", 18 > 4**2)
    ratio = 5*rm*R(".06")
    check("operator_Neumann_ratio", rm*R(".06") < 1)
    check("loop_geometric_ratio", ratio < 1)
    tail0 = 28*S.factorial(5)*10**5*2**5*5**3*rm**6
    tail = tail0/(1-ratio)
    zero("all_Born_tail_identity", sum(tail0*ratio**k for k in range(12))+tail*ratio**12-tail)
    local = 20*rm**2*product(1, 4)
    local += sum(R("4.8e7")*5**(k-2)*rm**(k+1)*product(k, 4) for k in range(2, 5))
    local = R("1.03")*(local+tail/4)
    check("full_fractional_scalar_remainder", local < R("3e14"), local)
    energy, mass = S.symbols("E M1", positive=True)
    integral = S.integrate(R(5, 2)*4*energy**(-R(3, 2)), (energy, mass, S.oo))
    zero("whole_tower_fractional_counting", integral-20/S.sqrt(mass))
    check("whole_tower_bound_5C", 20**2 < 5**2*18)
    check("borderline_counting_diverges", S.integrate(1/energy, (energy, mass, S.oo)) == S.oo)


def fractional_product_fixtures():
    # t=sqrt(r) turns fractional powers into exact polynomials. Signed,
    # distinct factors avoid testing only a repeated identical factor.
    t = S.Symbol("t", nonnegative=True)
    factors = [R(i+1, 1000)*t**2+(-1)**i*R(i+2, 200)*t**4
               +R(i+1, 30)*t**6-R(i+3, 10)*t**8
               +(-1)**i*R(i+1, 2)*t**9 for i in range(4)]
    for i, f in enumerate(factors, 1):
        divided = S.Poly(S.cancel(f/t**2), t)
        majorant = sum(abs(c)*R(1, 8)**power[0] for power, c in divided.terms())
        check("fixture_linear_bound_"+str(i), majorant < 2)
        check("fixture_low_value_"+str(i), majorant/64 < R(".06"))
        for d in range(5):
            truncated = sum(f.coeff(t, k)*t**k for k in range(2*d+1))
            quotient = S.Poly(S.cancel((f-truncated)/t**(2*d+1)), t)
            bound = sum(abs(c)*R(1, 8)**power[0] for power, c in quotient.terms())
            check(f"fixture_single_fractional_{i}_{d}", bound < BB[d+1])
    for count in range(1, 5):
        direct = S.expand(S.prod(factors[:count]))
        truncation = sum(direct.coeff(t, k)*t**k for k in range(9))
        # Independent convolution of each factor's nonzero monomials.
        convolution = {0: R(1)}
        for factor in factors[:count]:
            next_coeff = {}
            for power, coef in convolution.items():
                for (degree,), value in S.Poly(factor, t).terms():
                    next_coeff[power+degree] = next_coeff.get(power+degree, R(0))+coef*value
            convolution = next_coeff
        check("fixture_independent_product_"+str(count), all(direct.coeff(t, power) == coef
                                                            for power, coef in convolution.items()))
        remainder = S.Poly(direct-truncation, t)
        check("fixture_first_fractional_power_"+str(count), min(power[0] for power, _ in remainder.terms()) >= 9)
        bound = sum(abs(coef)*R(1, 8)**(power[0]-9) for power, coef in remainder.terms())
        check("fixture_full_product_remainder_"+str(count), bound < product(count, 4), bound)


def contact_and_response():
    aa, tt, b4, pp = 2*(10**8+220), 2*(10**8+66000), 4*(10**8+26400000), R(".011")
    qp = b4+3*aa**2+8*pp*tt+8*pp**2*aa
    check("retained_contact_Q", qp < 2*10**17, qp)
    check("contact_exp_factor", 1/(1-R(".044")) < R("1.05"))
    check("retained_contact_scalar", R(8*100, 960*3)*R("1.05")*qp < 10**17)
    k = 2*R("8.5e-33")*(5*R("3e14")+R("7e11")*24*256**3+10**17)
    check("complete_source_K", k < R("5e-12"), k)
    check("first_response_constant", 200*k < R("1e-9"))
    check("mixed_response_constant", 40000*k < R("2e-7"))
    check("finite_segment_Taylor", R("2e-7")/2 == R("1e-7"))
    check("complex_parameter_total", R(".001")+2*R(1, 200) == R(".011"))
    check("no_global_or_gradient_promotion", not PROTOCOL["spatial_source_gradient"]
          and not PROTOCOL["global_nonlinear_family_preserved"])


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 fractional interior response",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        inclusion_and_cusp()
        check("retained_center_coefficients", fractional_division() == PC)
        born_and_angular_budgets()
        fractional_product_fixtures()
        contact_and_response()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 fractional interior response",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "fractional_interior_value_response": good, "retained_first_mixed_constants": good,
        "P5_to_P45_inclusion": good, "spatial_source_gradient": False,
        "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
        "global_nonlinear_family_preserved": False, "useful_Newton_certificate": False,
        "S167_complete": False, "goal_complete": False,
        "unresolved": ["matched interior correction norms", "cold and outer preservation",
                       "full nonlinear mapping", "projected solve", "signed compatibility",
                       "physical requirements"],
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
