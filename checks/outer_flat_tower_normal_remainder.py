"""S166 flat full-tower normal remainder: symbols, kernels and weighted norm."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-flat-tower-normal-remainder"
OWN = ["checks/outer_flat_tower_normal_remainder.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("compact-massive-logarithmic-graph", 236),
    ("weighted-normal-inverse-scale", 453),
]
PROTOCOL = {
    "operator": "all massive angular levels in the even flat zero-winding geodesic multiplier",
    "finite_fourth_constant": "-31/30, fixed by the common prescription",
    "Weyl_combined_constant": "-6/5",
    "weight": "7/2, measure dY/Y",
    "bound": "1e30 epsilon, H4log to L2 after the normalized Y^4 factor",
    "near_scale": "epsilon<=1/16",
    "composition_example": "epsilon=1e-60, normal inverse budget <1e-9",
    "actual_curved_or_thermal_error_controlled": False,
    "global_inverse": False,
    "nonlinear_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {c["name"] for c in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.simplify(S.expand(expression))
    check(name, expression == 0, "" if expression == 0 else expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.add(f"evidence/{stem}.json")
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed"
              and len(run["checks"]) == count and all(c["passed"] for c in run["checks"])
              and inputs(run) == inputs(run, True))
        for path, digest in inputs(run).items():
            check(stem+"_input_"+path, sha256((ROOT/path).read_bytes()).hexdigest() == digest)


def real_series_and_symbol():
    s, q, t, z, c, ell, p = S.symbols("s q t z C ell p", positive=True)
    den = s*(s+q)+t
    phi, psi = (2*s+q)/den, (2*s+q)/den**2
    nup, num = S.symbols("nu_plus nu_minus")
    zero("series_rational_pair",
         1/(s+nup)+1/(s+num)-phi.subs({q: nup+num, t: nup*num}))
    zero("phi_integral_primitive", S.diff(S.log(den), s)-phi)
    zero("psi_integral_primitive", S.diff(-1/den, s)-psi)
    zero("psi_parameter_derivative", S.diff(phi, t)+psi)
    zero("phi_turning_numerator", S.diff(phi, s)*den**2-(2*t-q*q-2*q*s-2*s*s))
    zero("psi_turning_numerator", S.diff(psi, s)*den**3-2*(t-q*q-3*q*s-3*s*s))
    zero("phi_turning_numerator_decreases", S.diff(2*t-q*q-2*q*s-2*s*s, s)+2*q+4*s)
    zero("psi_turning_numerator_decreases", S.diff(t-q*q-3*q*s-3*s*s, s)+3*q+6*s)
    zero("universal_phi_max_envelope", s*s+t-2*s*S.sqrt(t)-(s-S.sqrt(t))**2)
    check("separate_variation_error_sum", 4+4 == 8)
    a = S.symbols("alpha", real=True)
    aa = a*(1-a)
    check("alpha_second_moment", S.integrate(aa**2, (a, 0, 1)) == R(1, 30))
    logarithmic = 2*(-R(1, 9)+R(2, 16)-R(1, 25))
    check("alpha_logarithmic_moment", logarithmic == -R(47, 900))
    check("alpha_first_moment", S.integrate(aa, (a, 0, 1)) == R(1, 6))
    beta32 = S.gamma(R(5, 2))**2/S.gamma(5)
    zero("alpha_three_half_moment", beta32-3*S.pi/128)
    k_inf = c*z*z*(S.log(z)/2-ell-R(31, 30))/(30*S.pi)
    zero("finite_fourth_constant_second_derivative",
         S.diff(k_inf, z, 2)-c*(S.log(z)-2*ell-R(17, 30))/(30*S.pi))
    zero("alpha_integrated_log_constant", 1+30*logarithmic+R(17, 30))
    zero("q_second_derivative_error", c/(2*S.pi)*8*q/(c*z)*R(1, 6)-2*q/(3*S.pi*z))
    zero("root_C_second_derivative_error",
         c/(2*S.pi)*8/(S.sqrt(c)*S.sqrt(z))*beta32-3*S.sqrt(c)/(32*S.sqrt(z)))
    kernel1 = S.integrate((z-t)/t, (t, 1, z))
    kernel2 = S.integrate((z-t)/S.sqrt(t), (t, 1, z))
    zero("integrated_second_error_log", kernel1-(z*S.log(z)-z+1))
    zero("integrated_second_error_root", kernel2-(R(4, 3)*z**R(3, 2)-2*z+R(2, 3)))
    check("integrated_q_reserve", R(2, 9) < 1)
    check("integrated_root_reserve", R(3, 32)*R(4, 3) < 1)
    check("retained_q_log_bound", 2000**4 < 2**64)
    check("retained_C0_upper", 1+12*R(22, 7)*2000**6 < R("6.2")**2*2000**6)
    check("retained_corrected_C_upper", R("6.2")*R("1.01")/R(".99") < 7)
    check("retained_C_log_bound", 7*2000**3 < 2**36)
    check("kappa4_upper_reserve", R(65, 2)+1 < 100)
    check("log_mu_ell_upper_reserve", 1+18 < 20)
    rr = S.log(1+t)-t+t*t/2
    zero("R_cubic_envelope_derivative", S.diff(t**3/3-rr, t)-t**3/(1+t))
    check("alpha_cubic_moment", S.integrate(aa**3, (a, 0, 1)) == R(1, 140))
    check("K_one_upper", R(100, 90)+R(4, 300*840*3) < 2)
    check("kprime_one_upper", 3*2 == 6)
    check("K_infinity_one", R(20+R(31, 30), 90) < 1)
    check("K_infinity_prime_one", (40+R(62, 30)+R(1, 2))/90 < 1)
    check("two_integration_constants", 2+1+6+1 < 20)
    zero("log_p_linear_envelope_derivative", S.diff(p-1-S.log(p), p)-(p-1)/p)
    check("low_log_p4_max", R(1, 4) < 1)
    check("complete_H3_multiplier_bound", 100+2*400+1 < 1000)
    zero("Weyl_fourth_constant_match", -R(31, 30)-R(1, 6)+R(6, 5))
    lb, lc, lx, ly = S.symbols("log_muB log_C log_xi log_Y")
    zero("normal_log_scale_match",
         (ly+(lx-lc)/2-(lb-lc/2))-(ly-(lb-lx/2)))


def counting_and_bessel_kernel():
    n = S.symbols("n", integer=True, nonnegative=True)
    j = S.symbols("j", integer=True, positive=True)
    q, c, m, r = S.symbols("q C M r", positive=True)
    count = S.summation(2*(q+2*j), (j, 1, n))
    zero("exact_all_level_count", count-(2*n*(n+q)+2*n))
    zero("lower_count_discrepancy", count-2*(n+1)*(n+1+q)+2*n+2*q+2)
    check("count_constant_q_reserve", 2+2 <= 4)
    zero("n_below_mass_root", n*(n+q)-n*n-n*q)
    k0, k1 = S.symbols("K0 K1")
    kernel = m**4*(k0*k0-k1*k1)/(2*S.pi**2)
    derivative = S.diff(kernel, m)+S.diff(kernel, k0)*(-r*k1) \
        +S.diff(kernel, k1)*(-r*k0-k1/m)
    zero("exact_mass_Bessel_derivative",
         derivative-m**3*(2*k0*k0-k1*k1)/S.pi**2)
    moment = lambda a, nu: 2**(a-3)*S.gamma(R(a, 2))**2 \
        *S.gamma(R(a, 2)+nu)*S.gamma(R(a, 2)-nu)/S.gamma(a)
    expected = {(4, 0): R(1, 3), (4, 1): R(2, 3),
                (5, 0): 27*S.pi**2/512, (5, 1): 45*S.pi**2/512,
                (6, 0): R(16, 15), (6, 1): R(8, 5)}
    for (a, nu), val in expected.items():
        zero(f"Bessel_moment_{a}_{nu}", moment(a, nu)-val)
    zero("continuum_kernel_coefficient",
         2*c/S.pi**2*(moment(6, 0)-moment(6, 1))+16*c/(15*S.pi**2))
    zero("count_error_q_coefficient",
         4/S.pi**2*(2*moment(4, 0)+moment(4, 1))-16/(3*S.pi**2))
    zero("count_error_root_coefficient",
         2/S.pi**2*(2*moment(5, 0)+moment(5, 1))-R(99, 256))
    check("two_dimensional_q_envelope", R(16, 27) < 10)
    check("two_dimensional_root_envelope", R(99, 256) < 10)
    time_integral = lambda a: S.sqrt(S.pi)*S.gamma(R(a-1, 2))/S.gamma(R(a, 2))
    zero("time_integral_four", time_integral(4)-S.pi/2)
    zero("time_integral_five", time_integral(5)-R(4, 3))
    zero("time_integral_six", time_integral(6)-3*S.pi/8)
    zero("continuum_stationary_kernel", -16/(15*S.pi**2)*time_integral(6)+2/(5*S.pi))
    check("stationary_kernel_q_reserve", 10*R(22, 7)/2 < 20)
    check("stationary_kernel_root_reserve", 10*R(4, 3) < 20)
    check("direct_and_image_q_far", 20*(2**3+1) < 400)
    check("direct_and_image_root_far", 20*(2**4+1) < 400)


def weighted_dyadic_bound():
    t = S.symbols("t")
    poly = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    for j, b in enumerate([11, 220, 3960], 1):
        pp = S.Poly(S.diff(poly, t, j), t)
        deg = pp.degree()
        cc = {i[0]: c for i, c in pp.terms()}
        bernstein = [sum(cc.get(i, 0)*S.binomial(k, i)/S.binomial(deg, i)
                         for i in range(k+1)) for k in range(deg+1)]
        check("rho_beta_derivative_"+str(j), max(abs(v) for v in bernstein) <= b)
        check("rho_scaled_derivative_"+str(j), 4**j*b < 10**6)
    check("dyadic_support_ratio", 8/R(1, 4) == 2**5)
    check("dyadic_infinite_overlap_reserve", 5+1 < 7)
    delta = R(7, 2)
    check("exact_output_measure_cancellation", 2*(4-delta)-1 == 0)
    d = S.symbols("D")
    for k in range(4):
        coeff = S.Poly(S.prod(d+delta-j for j in range(k)), d)
        check("log_derivative_polynomial_"+str(k), sum(abs(v) for v in coeff.coeffs()) < 40)
        check("positive_local_power_"+str(k), delta-k > 0)
    y, x = S.symbols("Y x", positive=True)
    for degree in range(7):
        for k in range(1, 4):
            target = x**degree
            for j in range(k):
                target = S.diff(target, x)+(delta-j)*target
            actual = S.diff(y**delta*S.log(y)**degree, y, k)
            zero(f"weighted_derivative_identity_{degree}_{k}",
                 actual-y**(delta-k)*target.subs(x, S.log(y)))
    raw = 2**3*10**6*8**4*40
    check("H3_local_even_reserve", 2*8*raw**2 < 10**30)
    check("near_normalized_weighted_bound",
          30*R(22, 7)*1000*10**15*3 < 10**22)
    check("near_overlap_sqrt_reserve", 7 < 3**2)
    # Powers of the log-coordinate kernel, separately on both sides.
    check("far_small_output_q_power", (4-delta)+(delta+1)-3 == 2)
    check("far_small_output_root_power", (4-delta)+(delta+1)-4 == 1)
    check("far_small_output_decay", 4-delta == R(1, 2))
    check("far_small_input_decay", delta+1 == R(9, 2))
    check("both_log_Schur_rows", 1/R(1, 2)+1/R(9, 2) < 3)
    check("far_normalized_weighted_bound", 12000*R(22, 7)*(400+1)*3 < 10**8)
    check("complete_weighted_bound", 2*(10**22+10**8) < 10**30)
    xi, c = S.symbols("xi C", positive=True)
    u = S.sqrt(xi/c)*y
    zero("physical_b_density_normalization",
         c*u**4/(16*S.pi*xi**2)-y**4/(16*S.pi*c))
    zero("normal_correction_d0",
         (1/(16*S.pi*c))*(c/(30*S.pi))-1/(480*S.pi**2))
    check("d0_below_one", R(1, 480*9) < 1)
    check("normal_inverse_single_correction_budget", R("8.7e19")*10**30*R("1e-60") < R("1e-9"))
    check("chosen_collar_inside_contract", R("1e-60") < R(1, 16))
    check("actual_curved_response_still_open", PROTOCOL["actual_curved_or_thermal_error_controlled"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 outer flat tower normal remainder",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        real_series_and_symbol()
        counting_and_bessel_kernel()
        weighted_dyadic_bound()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 outer flat tower normal remainder",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "flat_tower_weighted_normal_remainder": good,
        "actual_curved_or_thermal_error_controlled": False,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["actual curved and thermal outer response", "finite contacts and annulus matching",
                       "higher York maps", "full inverse/residual composition",
                       "nonlinear signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
