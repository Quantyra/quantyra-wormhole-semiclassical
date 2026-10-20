"""S157: shrinking charts and fixed-chart regular outer response."""
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
STEM = "comparison-path-outer-regular-response"
INPUTS = [
    "checks/comparison_path_outer_regular_response.py",
    "notes/comparison-path-outer-regular-response.md",
    "notes/comparison-path-interior-response.md",
    "evidence/comparison-path-interior-response.json",
    "notes/comparison-path-cold-response.md",
    "evidence/comparison-path-cold-response.json",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-centered-conormal-source.md",
    "notes/divided-jet-dirac-response.md",
    "notes/coincident-normal-dirac-response.md",
    "notes/outer-ward-and-constraint-transfer.md",
]
PROTOCOL = {
    "digits": 70, "scalar_y": [".08", ".01", "exp(-20)"],
    "parameters": ["0", ".4", "1"], "contour_nodes": 128,
    "contour_radius": ".005", "contour_tolerance": "1e-25",
    "scalar_steps": ["1e-7", "5e-8"], "scalar_difference_tolerance": "1e-6",
    "pressure_y": [".001", ".02", ".06"], "gauss_nodes": [32, 64],
    "pressure_steps": [".001", ".0005"], "pressure_normalization": "1e-14",
    "pressure_difference_tolerance": "1e-7", "quadrature_tolerance": "1e-24",
    "refinement_ratio": ".35", "refinement_floor": "1e-35",
    "endpoint_differences": "second-order one-sided within [0,1]",
    "fixtures_are_quantum_sources": False,
}
ROWS, NUMERICS = [], []


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
    X, r = s.symbols("X r", positive=True)
    yy = s.Rational(2, 25)
    H = 20*X/(1-20*yy*X)
    remainder = s.series(sum(H**k for k in range(6)), X, 0, 6).removeO().expand()
    cubic = s.series((s.Rational(21, 20)+H/100)**3, X, 0, 6).removeO().expand()
    Lmax = mp.mpf(1)-mp.log(mp.mpf(".08"))
    totals = []
    for j in range(2, 6):
        val = (mp.mpf("4e-8")*mp.mpf(".08")**4*Lmax*mp.mpf(str(remainder.coeff(X, j)))
               +mp.mpf("1e-11")*mp.mpf(".08")**3*mp.mpf(str(cubic.coeff(X, j))))
        check(f"scaled_metric_jet_{j}", val < mp.mpf("2.6e-5"), mp.nstr(val, 30))
        totals.append(mp.nstr(val, 30))
    for j in range(1, 6):
        check(f"seed_coordinate_cauchy_coefficient_{j}",
              Q(".01")*(j+1)*1000**j <= 20*2000**(j-1))
    check("seed_complex_coordinate_margin", Q(".0042")*Q("1.12") < Q(".0048") < Q(".01"))
    check("real_cap_proper_margin", Q("1.03")/40 < Q(".05"))
    check("scaled_norm_below_existing_reserve", Q("2.6e-5") < Q("1.3e-4"))
    check("ellipse_stays_in_scaled_cone", Q("64.5009765625")*Q("1.3e-4") < Q(".011"))
    check("proper_cubic_coefficient", Q("4e-11")*Q(".075") < Q("1e-11"))
    check("small_mass_evaluation_fits", Q(1, 320) < Q(1, 160))
    M = s.symbols("M", positive=True)
    check("two_row_small_mass_factor",
          s.expand(M*(1+32/(M*r))**2-(M+64/r+1024/(M*r**2))) == 0)
    check("row_commutator_scalar_constant", 2*32000**2*320*Q("1.02") < 9*10**11)
    # Homogeneity follows through the full finite product recursion.
    p = [None, s.Rational(2, 100), s.Integer(300), s.Integer(50000), s.Integer(7000000)]
    B = [None, 2, 400, 50000, 7000000, 2000000000]
    def product(n, d, scale):
        if n == 1:
            return B[d+1]/scale**(d+1)
        if d == n-1:
            return (2/scale)**n
        return (B[d-n+2]/scale**(d-n+2)*(2/scale)**(n-1)
                +sum(p[a]/scale**a*product(n-1, d-a, scale) for a in range(1, d-n+2)))
    for n in range(1, 5):
        check(f"born_degree_five_scaling_{n}", s.simplify(product(n, 4, r)*r**5-product(n, 4, s.Integer(1))) == 0)
    check("born_neumann_ratio", 5*Q("1.02")*Q(".06") < 1)
    eta = Q("8.5e-33")
    angular = 24*320**3+64*8*320**2+1024*8*320
    Kunit = 2*eta*(Q(2, 3)*10**17+10**17+9*10**11*angular)
    check("complete_shrinking_complex_bound", Kunit < Q("2e-11"), str(Kunit))
    check("difference_and_cubic_complex_budget", 2*Q("2e-11")+Q("1e-30") < Q("1e-10"))
    E, delta = s.symbols("E delta", positive=True)
    check("angular_count_zero_moment", s.integrate(4*delta*E**2*s.exp(-delta*E), (E, 0, s.oo)) == 8/delta**2)
    check("angular_count_inverse_moment", s.integrate(4*(1+delta*E)*s.exp(-delta*E), (E, 0, s.oo)) == 8/delta)
    a = Q("1.2e-8")
    H0, H1, H2 = Q("1.01"), Q("5e-8"), Q("2.4e-15")
    L0, L1, L2 = Q(".01"), Q("8e-10"), Q("5.5e-17")
    I = Q(".08")**2/2
    check("ward_H_first", 4*a*H0 < H1)
    check("ward_H_second", 16*a*a*H0 < H2)
    check("ward_L_value", H0*3*I < L0)
    check("ward_L_first", H0*(12*a+Q("1e-7"))*I < L1)
    check("ward_L_second", H0*(48*a*a+8*a*Q("1e-7"))*I < L2)
    lowerL = Q("3.5")
    P0, P1, P2 = Q("6e-10"), Q("5.2e-13"), Q("1.1e-14")
    S0, S1, S2 = Q("3e-7"), Q("6e-14"), Q("9e-14")
    bound0 = H0*P0/lowerL**2+L0*S0
    bound1 = (H0*P1+H1*P0)/lowerL**4+L0*S1+L1*S0/lowerL**2
    bound2 = (H0*P2+2*H1*P1+H2*P0)/lowerL**6+L0*S2+2*L1*S1/lowerL**2+L2*S0/lowerL**4
    for j, v, bound in zip(range(3), [bound0, bound1, bound2], [Q("4e-9"), Q("5e-15"), Q("1e-15")]):
        check(f"pressure_weighted_budget_{j}", v < bound, str(v))
    N0, N1, N2 = Q("3.1e-7"), Q("7e-14"), Q("9.2e-14")
    check("regular_profile_value", 2*Q("4e-9")+S0 < N0)
    check("regular_profile_first", 2*Q("5e-15")+S1 == N1)
    check("regular_profile_second", 2*Q("1e-15")+S2 == N2)
    dg = Q("3e-7")
    m2 = Q("1.003")*N2+Q("4.02")*dg*N1+Q("5.54")*dg**2*N0
    a2 = Q("1.006")*N2+Q("10.08")*dg*N1+Q("26.5")*dg**2*N0
    outs = [Q(".302")*m2/2, 1580*m2/2, Q(".019")*a2/2]
    for tag, v, bound in zip(["mass", "density", "lapse"], outs, [Q("1.4e-14"), Q("7.3e-11"), Q("9e-16")]):
        check("outer_regular_"+tag, v < bound, str(v))
    check("outer_cross_terms_retained", m2 > Q("1.003")*N2 and a2 > Q("1.006")*N2)
    NUMERICS.append({"scaled_jets": totals, "constraint_second_bounds_exact": [str(x) for x in outs],
                     "full_outer_source_bounded": False, "goal_complete": False})


def numerical_budgets():
    y = mp.mpf(".08")
    L = 1-mp.log(y)
    check("seed_strip_log_guard", mp.mpf("1.01")*mp.exp(mp.mpf(".09")) < mp.mpf("1.12"))
    check("capped_gradient", mp.exp(mp.mpf(".027"))*(mp.mpf("8e-6")+160*mp.mpf("3e-8")*y*y) < mp.mpf("1e-5"))
    check("radial_direction_low_guard", 3*mp.mpf("1e-11")*y+mp.mpf("2e-6")*y*y*L < mp.mpf("1e-7"))
    check("real_collar_log_lower", L > mp.mpf("3.5"))
    check("outer_real_endpoint_budget", 2*mp.mpf("8.5e-33")*mp.mpf("1e24")*(mp.mpf("1e-7")+mp.mpf("4e-11"))*mp.mpf("1.01") < mp.mpf("5e-15"))
    check("compact_real_endpoint_budget", mp.mpf("3.4e-14")+mp.mpf("1e-20") < mp.mpf("5e-15")*L**2)
    check("seed_regular_value_budget", 2*mp.mpf("8.5e-33")*400**2*mp.mpf("1e20")+mp.mpf("5e-15") < mp.mpf("3e-7"))
    check("real_to_complex_log_bound", mp.log(20000) < 7+2*mp.log(1/y))
    check("outer_scalar_first_budget", mp.mpf("5e-15")*mp.exp(2)*(7/mp.log(256))**2 < mp.mpf("6e-14"))
    check("outer_scalar_second_budget", mp.mpf("5e-15")*mp.exp(4)/8*(7/mp.log(256))**4 < mp.mpf("9e-14"))
    check("outer_region_covered", mp.sqrt(8)*mp.log(mp.mpf(40)/39) < y)
    check("outer_log_weight_conversion", mp.sqrt(mp.mpf("5.9"))/mp.mpf("1.001") > 1)
    ell = 1+mp.log(40)
    I6 = sum(mp.binomial(6, j)*mp.factorial(j)*ell**(6-j) for j in range(7))/40
    J6 = sum(mp.binomial(6, j)*mp.factorial(j)*ell**(6-j)/4**(j+1) for j in range(7))/40**4
    check("sixth_log_mass_integral", I6 < 1580, mp.nstr(I6, 30))
    check("sixth_log_primitive_constant", mp.mpf("12.2")*I6/40**3 < mp.mpf(".302"))
    check("sixth_log_lapse_constant", mp.mpf("12.2")*J6 < mp.mpf(".019"), mp.nstr(mp.mpf("12.2")*J6, 30))
    for power, exact in [(1, I6), (4, J6)]:
        value = mp.quad(lambda z: (ell+z)**6*mp.exp(-power*z), [0, 1, 4, 16, mp.inf])/40**power
        check(f"sixth_log_independent_integral_{power}", abs(value-exact) < mp.mpf("1e-55"))
    NUMERICS.append({"I6": mp.nstr(I6, 35), "J6": mp.nstr(J6, 35)})


def fd(fun, x, h, j):
    if x == 0 or x == 1:
        h = h if x == 0 else -h
        assert all(0 <= x+n*h <= 1 for n in range(4))
        v = [fun(x+n*h) for n in range(4)]
        return (-3*v[0]+4*v[1]-v[2])/(2*h) if j == 1 else (2*v[0]-5*v[1]+4*v[2]-v[3])/h**2
    return (fun(x+h)-fun(x-h))/(2*h) if j == 1 else (fun(x+h)-2*fun(x)+fun(x-h))/h**2


def differences(tag, fun, x, target, j, steps, tol):
    errors = [abs(fd(fun, x, mp.mpf(h), j)-target) for h in steps]
    check(tag+"_difference", max(errors) < mp.mpf(tol)*(1+abs(target)),
          {"target": mp.nstr(target, 25), "errors": [mp.nstr(e, 12) for e in errors]})
    check(tag+"_refinement", errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-35"))


def scalar_fixtures():
    x = s.symbols("x")
    for label, y in zip(PROTOCOL["scalar_y"], [mp.mpf(".08"), mp.mpf(".01"), mp.exp(-20)]):
        L = 1-mp.log(y)
        a, B = mp.mpf("5e-15")*L**2, mp.mpf("1e-10")*y**-5
        degree = int(mp.floor(mp.log(B/a)/mp.log(256)))
        check("fixture_ellipse_bound_"+label,
              (a/4)*(mp.mpf(256)**degree+mp.mpf(256)**(-degree)) <= B)
        poly = s.chebyshevt(degree, 2*x-1)/2
        derivs = [s.lambdify(x, s.diff(poly, x, j), "mpmath") for j in [1, 2]]
        fun = lambda q: mp.chebyt(degree, 2*q-1)/2
        for parameter in PROTOCOL["parameters"]:
            lam = mp.mpf(parameter)
            for j in [1, 2]:
                tag = f"scalar_{label}_{parameter}_{j}"
                target = derivs[j-1](lam)
                bound = mp.mpf("6e-14")*L**4/a if j == 1 else mp.mpf("9e-14")*L**6/a
                check(tag+"_weighted_bound", abs(target) < bound)
                radius, count = mp.mpf(".005"), 128
                contour = mp.factorial(j)*sum(
                    fun(lam+radius*mp.exp(2j*mp.pi*n/count))*mp.exp(-2j*mp.pi*j*n/count)
                    for n in range(count))/count/radius**j
                check(tag+"_contour", abs(contour-target) < mp.mpf("1e-25")*(1+abs(target)))
                differences(tag, fun, lam, target, j, PROTOCOL["scalar_steps"], "1e-6")
        NUMERICS.append({"scalar_y": label, "chebyshev_degree": degree, "amplitude_normalized": True})


def pressure_fixtures():
    nodes = {n: mp.gauss_quadrature(n, "legendre") for n in [32, 64]}
    top, normalization = mp.mpf(".08"), mp.mpf("1e-14")
    def u(v):
        return mp.mpf("1e-8")*v**3+mp.mpf("2e-8")*v**4*(1-mp.log(v))
    def up(v):
        return mp.mpf("3e-8")*v**2+mp.mpf("2e-8")*v**3*(3-4*mp.log(v))
    def source(v, lam, j=0):
        L = 1-mp.log(v)
        coeff = [mp.mpf("3e-10")+lam*mp.mpf("2e-13")+lam**2*mp.mpf("4e-15"),
                 mp.mpf("2e-13")+lam*mp.mpf("8e-15"), mp.mpf("8e-15")][j]
        return coeff*L*L
    def memory(lam, j=0):
        return [mp.mpf("5e-10")+lam*mp.mpf("4e-13")+lam**2*mp.mpf("5e-15"),
                mp.mpf("4e-13")+lam*mp.mpf("1e-14"), mp.mpf("1e-14")][j]
    for label in PROTOCOL["pressure_y"]:
        y = mp.mpf(label)
        length = mp.log(top/y)
        def integrate(fun, count):
            xs, ws = nodes[count]
            return sum(w*length/2*fun(y*mp.exp(length*(z+1)/2))*y*mp.exp(length*(z+1)/2)
                       for z, w in zip(xs, ws))
        def pressure(lam, count):
            H = mp.exp(top**2-y**2+2*lam*(u(top)-u(y)))
            return H*memory(lam)-integrate(lambda v:
                mp.exp(v*v-y*y+2*lam*(u(v)-u(y)))*(v+lam*up(v))*source(v, lam), count)
        for parameter in PROTOCOL["parameters"]:
            lam = mp.mpf(parameter)
            dt = u(top)-u(y)
            H = mp.exp(top**2-y**2+2*lam*dt)
            mem2 = H*(memory(lam, 2)+4*dt*memory(lam, 1)+4*dt*dt*memory(lam))
            def density(v, omit=False):
                d = u(v)-u(y)
                K = mp.exp(v*v-y*y+2*lam*d)
                A = v+lam*up(v)
                ans = A*source(v, lam, 2)
                if not omit:
                    ans += (4*d*A+2*up(v))*source(v, lam, 1)+(4*d*d*A+4*d*up(v))*source(v, lam)
                return K*ans
            integral = mp.quad(lambda z: density(y*mp.exp(z))*y*mp.exp(z), [0, length/2, length])
            target = (mem2-integral)/normalization
            tag = f"pressure_{label}_{parameter}"
            for count in [32, 64]:
                check(tag+f"_second_quadrature_{count}",
                      abs((mem2-integrate(density, count))/normalization-target) < mp.mpf("1e-24"))
            check(tag+"_node_agreement", abs(pressure(lam, 32)-pressure(lam, 64))/normalization < mp.mpf("1e-24"))
            differences(tag, lambda q: pressure(q, 64)/normalization, lam, target, 2,
                        PROTOCOL["pressure_steps"], "1e-7")
            omitted = (H*memory(lam, 2)-integrate(lambda v: density(v, True), 64))/normalization
            check(tag+"_kernel_omission_detected", abs(target-omitted) > mp.mpf("1e-13"))
            check(tag+"_memory_omission_detected", abs(mem2/normalization) > mp.mpf(".5"))
            NUMERICS.append({"pressure_y": label, "parameter": parameter,
                             "normalized_second": mp.nstr(target, 25)})


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = 70
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S157", "input_count": len(before)}), flush=True)
    try:
        predecessors()
        exact_checks()
        numerical_budgets()
        scalar_fixtures()
        pressure_fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S157", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error, "scope": "fixed-chart actual-path outer regular canonical response",
        "unresolved": ["allocated cubic/Weyl and other sources", "general family and full inverse",
                       "coupled solution", "physical completion requirements"],
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
