"""S155: fixed-chart transfer and actual comparison-path response."""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "comparison-path-interior-response"
INPUTS = [
    "checks/comparison_path_interior_response.py",
    "notes/comparison-path-interior-response.md",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/divided-jet-dirac-response.md",
    "evidence/divided-jet-dirac-response.json",
    "notes/interior-pressure-and-bolt-primitive.md",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/coincident-normal-dirac-response.md",
    "notes/curved-remote-dirac-response.md",
]
PROTOCOL = {
    "digits": 70, "coordinate_frequencies": [200, 1000, 4000],
    "coordinate_points": [["0", "0"], [".01", "-.006"], [".005", "-.003"]],
    "coordinate_jet_tolerance": "1e-20", "coordinate_jet_normalization": "j!100^j",
    "pressure_backgrounds": [".3", ".8"], "pressure_cuts": [".1", "1", "4"],
    "pressure_steps": [".001", ".0005"], "gauss_nodes": [32, 64],
    "pressure_panel_max_length": ".25", "pressure_second_normalization": "1e-14",
    "pressure_difference_tolerance": "1e-7", "refinement_ratio": ".35",
    "refinement_floor": "1e-40", "pressure_quadrature_tolerance": "1e-25",
    "fixtures_are_quantum_sources": False,
}
ROWS, NUMERICS = [], []


def check(name, passed, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessor_checks():
    for stem, count in [("divided-jet-dirac-response", 93),
                        ("fixed-chart-global-constraint-transfer", 99)]:
        run = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))["runs"][-1]
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(r["passed"] for r in run["checks"]))
        check(stem+"_current_inputs", run["unchanged_inputs"] and all(
            sha256((ROOT/p).read_bytes()).hexdigest() == h
            for p, h in run["input_sha256_after"].items()))


def exact_checks():
    X, a, c = s.symbols("X a c", nonzero=True)
    inverse = sum((-c)**(j-1)*s.catalan(j-1)/a**(2*j-1)*X**j for j in range(1, 6))
    check("inverse_quadratic_coefficients",
          s.series(a*inverse+c*inverse**2-X, X, 0, 6).removeO().expand() == 0)
    check("inverse_catalan_first_terms", [s.catalan(j) for j in range(5)] == [1, 1, 2, 5, 14])
    amin = s.Rational(999999, 1000000)
    cc = s.Rational(4, 10**6)
    H = (1+s.Rational(1, 10**6))*X
    H += sum(2*s.catalan(j-1)*cc**(j-1)/amin**(2*j-1)*X**j for j in range(2, 6))
    D = s.Rational(5, 10**9)+H-X
    d = s.Rational(1, 160)
    composition = s.Poly(s.series(sum((100*H)**k for k in range(6)), X, 0, 6).removeO(), X)
    seed = s.Poly(s.series(s.Rational(27, 1000)*(
        (1-(X+D)/d)**-2-(1-X/d)**-2), X, 0, 6).removeO(), X)
    jacobian = s.Poly(s.series(sum(2*(2*cc*H/amin)**k/k for k in range(1, 6)),
                                   X, 0, 6).removeO(), X)
    displayed = [s.Rational(101086, 10**9), s.Rational(103131, 10**9),
                 s.Rational(107893, 10**9), s.Rational(118290, 10**9)]
    totals = []
    for j, upper in zip(range(2, 6), displayed):
        total = (s.Rational(1, 10000)*composition.coeff_monomial(X**j)
                 + seed.coeff_monomial(X**j)+jacobian.coeff_monomial(X**j))/100**j
        check(f"fixed_chart_divided_jet_{j}", total < upper < s.Rational(13, 10**5), str(total))
        totals.append(str(total))
    check("chart_image_margin", (Q(1)-Q(2, 10**8))/32-Q(4, 10**6)/32**2 > Q(1, 40))
    check("inverse_denominator", Q(1)-Q(2, 10**8)-2*Q(4, 10**6)/32 > Q(999999, 10**6))
    check("inverse_displacement", Q(2, 10**8)/32+Q(4, 10**6)/32**2 < Q(5, 10**9))
    check("inverse_linear_coordinate_norm",
          2*(Q(2, 10**8)+2*Q(4, 10**6)/32)/Q(999999, 10**6) < Q(1, 10**6))
    check("proper_middle_cover", Q(59, 10)/(40*Q(1001, 1000))**2 > Q(606, 10000)**2)
    check("point_evaluation_disk_fits", Q(1, 180) < Q(1, 160))
    pp = [None, Q(2, 100), Q(300), Q(50000), Q(7000000)]
    remainders = [None, 400, 50000, 7000000, 2000000000]
    for degree in range(1, 5):
        value = Q(202, 100)*80**degree + sum(pp[j]*80**(degree+1-j)
                                                            for j in range(2, degree+1))
        check(f"smaller_cap_taylor_division_{degree}", value < remainders[degree], str(value))
    check("actual_cap_rows_and_commutator",
          2*32000**2*320*Q(102, 100) < 9*10**11)
    eta = Q(85, 10**34)
    K = 2*eta*(Q(2, 3)*10**17+9*10**11*24*320**3+10**17)
    check("smaller_cap_complete_complex_bound", K < Q(125, 10**13), str(K))
    check("smaller_cap_first_constant", 200*K < Q(25, 10**10))
    check("smaller_cap_mixed_constant", 40000*K < Q(5, 10**7))
    U = Q(13, 10**5)
    S0, S1, S2 = Q(5, 10**10), Q(325, 10**15), Q(9, 10**15)
    check("path_background_inside_source_ball", U < Q(1, 1000))
    check("path_scalar_first", Q(25, 10**10)*U == S1)
    check("path_scalar_second", Q(5, 10**7)*U*U < S2)
    check("path_scalar_value", Q(4352, 10**13)+S1 < S0)
    low, grad = Q(12, 10**9), Q(3, 10**6)
    F1 = 2*Q(24, 10)*low+84*grad
    F2 = 8*Q(24, 10)*low**2+8*84*low*grad
    check("path_ward_first_kernel", F1 < Q(253, 10**6))
    check("path_ward_second_kernel", F2 < Q(25, 10**12))
    P1, P2 = Q(52, 10**14), Q(11, 10**15)
    check("path_pressure_first", Q(12, 10)*S1+Q(253, 10**6)*S0 < P1)
    check("path_pressure_second", Q(12, 10)*S2+2*Q(253, 10**6)*S1+Q(25, 10**12)*S0 < P2)
    N0, N1, N2 = Q(17, 10**10), Q(14, 10**13), Q(31, 10**15)
    check("path_profiles_value", 2*Q(12, 10)*S0+S0 == N0)
    check("path_profiles_first", 2*P1+S1 < N1)
    check("path_profiles_second", 2*P2+S2 == N2)
    dg = Q(3, 10**7)
    m2 = Q(1003, 1000)*N2+Q(402, 100)*dg*N1+Q(554, 100)*dg**2*N0
    a2 = Q(1006, 1000)*N2+Q(1008, 100)*dg*N1+Q(265, 10)*dg**2*N0
    check("middle_mass_second", 279*m2/2 < Q(44, 10**13))
    check("middle_infinity_second", m2/2 < Q(16, 10**15))
    check("middle_lapse_second", 1300*a2/2 < Q(21, 10**12))
    check("geometry_cross_terms_retained", m2 > Q(1003, 1000)*N2 and a2 > Q(1006, 1000)*N2)
    NUMERICS.append({"divided_jet_totals_exact": totals, "complex_K_exact": str(K),
                     "middle_second_bounds_exact": {"mass": str(279*m2/2),
                                                   "infinity": str(m2/2), "lapse": str(1300*a2/2)},
                     "full_geometry_is_solution": False, "general_family_preservation_proved": False})


def coordinate_checks():
    mp.mp.dps = PROTOCOL["digits"]
    check("normalizer_a_bound", mp.expm1(mp.mpf("1.2e-8")) < mp.mpf("2e-8"))
    check("normalizer_c_bound", mp.exp(mp.mpf("2.4e-8"))*mp.mpf("7.7e-6")/2 < mp.mpf("4e-6"))
    check("seed_log_bound", -mp.log(1-mp.mpf(".026")) < mp.mpf(".027"))
    check("proper_endpoint_cover", mp.exp(-mp.mpf("1.2e-8"))*mp.mpf(".051") > mp.mpf(".05"))
    check("local_low_gradient", mp.exp(mp.mpf(".027"))*mp.mpf("7.7e-6") < mp.mpf("8e-6"))
    check("capped_low_gradient", mp.exp(mp.mpf(".027"))*(
        mp.mpf("8e-6")+160*mp.mpf("1.2e-8")) < mp.mpf("1.1e-5"))
    check("smaller_cap_linear_potential", mp.exp(mp.mpf(".05"))*(
        8/(40*mp.mpf(".974"))+mp.mpf(".016")+mp.mpf(".011")*80) < 2)
    check("smaller_gluing_cutoff_slope", 310*mp.exp(mp.mpf(".027")) < 320)
    check("smaller_gluing_proper_separation", mp.exp(-mp.mpf(".05"))/160 > mp.mpf(1)/320)
    for freq in PROTOCOL["coordinate_frequencies"]:
        def b(x, y):
            return (mp.mpf("2e-10")+mp.mpf("1e-7")*x-mp.mpf("2e-7")*y
                    +mp.mpf("1e-10")*(mp.sin(freq*x)+mp.cos(freq*y))
                    +mp.mpf(".001")*abs(x-mp.mpf(".005"))**mp.mpf("5.5"))

        def sigma(x, y):
            return -mp.log1p(-(x*x+y*y)/4)

        a = mp.exp(-b(0, 0))
        bz = (mp.diff(lambda x: b(x, 0), 0)-mp.j*mp.diff(lambda y: b(0, y), 0))/2
        c = -a*a*bz
        check(f"complex_normalizer_{freq}", abs(mp.re(c)) > 0 and abs(mp.im(c)) > 0)

        def F(z):
            return a*z+c*z*z

        def G(z):
            return 2*z/(a+mp.sqrt(a*a+4*c*z))

        def normal(x, y):
            z = mp.mpc(x, y)
            f = F(z)
            return sigma(mp.re(f), mp.im(f))+b(mp.re(f), mp.im(f))+mp.log(abs(a+2*c*z))-sigma(x, y)

        def recovered(x, y):
            z = G(mp.mpc(x, y))
            return normal(mp.re(z), mp.im(z))+sigma(mp.re(z), mp.im(z))-sigma(x, y)-mp.log(abs(a+2*c*z))

        for index, pair in enumerate(PROTOCOL["coordinate_points"]):
            x, y = map(mp.mpf, pair)
            z = mp.mpc(x, y)
            root = mp.findroot(lambda t: F(t)-z, z/a)
            check(f"independent_inverse_{freq}_{index}", abs(root-G(z)) < mp.mpf("1e-50"))
            check(f"coordinate_value_identity_{freq}_{index}", abs(recovered(x, y)-b(x, y)) < mp.mpf("1e-50"))
            errors, jets = [], []
            for order in range(1, 6):
                for first in range(order+1):
                    multi = (first, order-first)
                    direct = mp.diff(b, (x, y), multi)
                    transformed = mp.diff(recovered, (x, y), multi)
                    scale = mp.factorial(order)*100**order
                    errors.append(abs(direct-transformed)/scale)
                    if order >= 2:
                        jets.append(abs(direct)/scale)
            check(f"coordinate_jets_through_five_{freq}_{index}",
                  max(errors) < mp.mpf("1e-20"), mp.nstr(max(errors), 12))
            check(f"fixture_high_jet_budget_{freq}_{index}", max(jets) < mp.mpf("1.3e-4"))
        omission = mp.log(abs(a+2*c*G(mp.mpc(".01", "-.006"))))
        check(f"jacobian_omission_detected_{freq}", abs(omission) > mp.mpf("1e-12"))


def pressure_checks():
    mp.mp.dps = PROTOCOL["digits"]
    nodes = {}
    for count in PROTOCOL["gauss_nodes"]:
        xs, ws = mp.gauss_quadrature(count, "legendre")
        nodes[count] = (list(xs), list(ws))

    def u(t):
        return mp.mpf("1e-8")*t*t/(1+t*t)

    def up(t):
        return mp.mpf("2e-8")*t/(1+t*t)**2

    def source(t, lam):
        return (mp.mpf("3e-10")+mp.mpf("1e-10")/(1+t*t)
                +lam*mp.mpf("3e-13")/(1+t*t)+lam*lam*mp.mpf("4e-15")/(1+t*t))

    def integrate(fun, cut, count):
        xs, ws = nodes[count]
        panels = int(mp.ceil(4*cut))
        total = mp.mpf(0)
        for j in range(panels):
            lo, hi = cut*j/panels, cut*(j+1)/panels
            mid, half = (lo+hi)/2, (hi-lo)/2
            total += half*sum(w*fun(mid+half*x) for x, w in zip(xs, ws))
        return total

    def pressure(lam, cut, count):
        denom = mp.sinh(cut)**2*mp.exp(2*lam*u(cut))
        def fun(t):
            w = mp.sinh(t)**2*mp.exp(2*lam*u(t))
            wp = 2*mp.sinh(t)*mp.cosh(t)*mp.exp(2*lam*u(t))+2*lam*up(t)*w
            return wp*source(t, lam)/(2*denom)
        return integrate(fun, cut, count)

    for background in PROTOCOL["pressure_backgrounds"]:
        lam = mp.mpf(background)
        for cuttext in PROTOCOL["pressure_cuts"]:
            cut = mp.mpf(cuttext)
            def second_density(t, omit=False):
                if t == 0:
                    return mp.mpf(0)
                d = u(t)-u(cut)
                kernel = (mp.sinh(t)/mp.sinh(cut))**2*mp.exp(2*lam*d)
                ell = 2*mp.coth(t)+2*lam*up(t)
                s1 = (mp.mpf("3e-13")+lam*mp.mpf("8e-15"))/(1+t*t)
                s2 = mp.mpf("8e-15")/(1+t*t)
                value = ell*s2
                if not omit:
                    value += (4*d*ell+4*up(t))*s1+(4*d*d*ell+8*d*up(t))*source(t, lam)
                return kernel*value/2
            exact = mp.quad(second_density, [0, cut/2, cut])/mp.mpf("1e-14")
            tag = background+"_"+cuttext
            for count in PROTOCOL["gauss_nodes"]:
                value = integrate(second_density, cut, count)/mp.mpf("1e-14")
                check(f"pressure_second_quadrature_{tag}_{count}",
                      abs(value-exact) < mp.mpf("1e-25"))
            check(f"pressure_value_node_agreement_{tag}",
                  abs(pressure(lam, cut, 32)-pressure(lam, cut, 64))/mp.mpf("1e-14") < mp.mpf("1e-25"))
            errors = []
            for steptext in PROTOCOL["pressure_steps"]:
                step = mp.mpf(steptext)
                value = (pressure(lam+step, cut, 64)-2*pressure(lam, cut, 64)
                         +pressure(lam-step, cut, 64))/step**2/mp.mpf("1e-14")
                errors.append(abs(value-exact))
            check(f"pressure_second_difference_{tag}",
                  max(errors) < mp.mpf("1e-7")*(1+abs(exact)),
                  {"normalized_target": mp.nstr(exact, 25), "errors": [mp.nstr(e, 12) for e in errors]})
            check(f"pressure_difference_refinement_{tag}",
                  errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-40"))
            omitted = mp.quad(lambda t: second_density(t, True), [0, cut/2, cut])/mp.mpf("1e-14")
            check(f"pressure_kernel_omission_detected_{tag}", abs(exact-omitted) > mp.mpf("1e-15"))
            NUMERICS.append({"pressure_background": background, "cut": cuttext,
                             "normalized_second": mp.nstr(exact, 25)})


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        predecessor_checks()
        exact_checks()
        coordinate_checks()
        pressure_checks()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S155", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error, "scope": "canonical interior and middle-constraint response on retained comparison path",
        "unresolved": ["general family admissibility", "cold/outer and other source mappings",
                       "full higher-order inverse and coupled solve", "physical completion requirements"],
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
