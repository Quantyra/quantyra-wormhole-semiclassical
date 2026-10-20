"""S166 compact-chart actual massive fourth/log graph verification."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "compact-massive-logarithmic-graph"
OWN = ["checks/compact_massive_logarithmic_graph.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-curved-quantum-tt", 157),
    ("strong-interior-quantum-remainder", 141),
    ("tt-mass-numerical-bound", 398),
]
PROTOCOL = {
    "operator": "actual canonical massive Hessian; optional bounded Weyl principal correction",
    "norm": "unitary R2 Fourier multiplier (1+|xi|^2)^2 log(2+|xi|)",
    "charts": "center Y>=1/4, input radius 1/256, observation radius 1/128",
    "arithmetic": "exact identities, rational continuum bounds, no fitted constants",
    "full_TT_inverse": False,
    "nonlinear_contraction": False,
    "low_norm_constant": "1e390, explicitly not a useful Newton certificate",
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {c["name"] for c in ROWS}
    item = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(item)
    print(json.dumps(item), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.factor(S.cancel(S.expand(S.simplify(expression))))
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
        out.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
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


def flat_multiplier():
    k = S.symbols("k", integer=True, positive=True)
    moment = 2**(2*k-2)*S.factorial(k)**4/(S.pi**2*k*S.factorial(2*k+1))
    angle = 2*S.pi*S.factorial(2*k)/(4**k*S.factorial(k)**2)
    flat = S.factorial(k)**2/(2*S.pi*k*S.factorial(2*k+1))
    zero("all_k_geodesic_flat_moment", moment*angle/S.factorial(2*k)-flat)
    # Independently derive the difference of the two squared-Bessel moments.
    j0 = 2**(2*k-1)*S.gamma(k+1)**4/S.gamma(2*k+2)
    j1 = 2**(2*k-1)*S.gamma(k+1)**2*S.gamma(k+2)*S.gamma(k)/S.gamma(2*k+2)
    zero("all_k_Bessel_square_difference", (j1-j0)/(2*S.pi**2)-moment)
    for n in range(3, 9):
        alpha = S.symbols("alpha")
        beta = S.integrate((alpha*(1-alpha))**n, (alpha, 0, 1))
        zero("independent_beta_moment_"+str(n), beta-S.factorial(n)**2/S.factorial(2*n+1))
    zero("sixth_flat_normalization", flat.subs(k, 3)-1/(840*S.pi))
    zero("fourth_restoration_normalization", 3*S.pi/(4*90*S.pi**2)-1/(120*S.pi))
    z = S.symbols("z", positive=True)
    rr = S.log(1+z)-z+z*z/2
    zero("R_positive_integrand", S.diff(rr, z)-z*z/(1+z))
    zero("R_scale_upper_derivative", S.diff(3*rr-z*S.diff(rr, z), z)-z**3/(1+z)**2)
    zero("R_scale_lower_derivative", S.diff(z*S.diff(rr, z)-2*rr, z)-z*z/(1+z)**2)
    zero("scale_upper_zero_origin", (3*rr-z*S.diff(rr, z)).subs(z, 0))
    zero("scale_lower_zero_origin", (z*S.diff(rr, z)-2*rr).subs(z, 0))
    aa, ss = S.symbols("a s", positive=True)
    zero("both_signs_scale_integration", 2*S.integrate(6*S.exp(6*ss), (ss, 0, aa))
         -2*(S.exp(6*aa)-1))
    # Proper normalization of the full directional sixth moment was checked above.
    M, p = S.symbols("M p", positive=True)
    zero("gap_moment_sum_order", (M**2*(p*p/M**2)**3)-p**6/M**4)


def exp_enclosure(value, n=24):
    low = sum(value**k/S.factorial(k) for k in range(n+1))
    tail = value**(n+1)/S.factorial(n+1)/(1-value/(n+2))
    return low, low+tail


def continuum_and_geometry():
    a = R(27, 1000)
    for name, value, upper in [("scale", 6*a, R(118, 100)),
                               ("Weyl", 4*a, R(112, 100)),
                               ("density", 2*a, R(11, 10))]:
        low, high = exp_enclosure(value)
        check("exponential_"+name, 1 < low < high < upper)
    low, _ = exp_enclosure(R(5, 4))
    check("log3_less_than_five_fourths", low > 3)
    low, _ = exp_enclosure(R(10, 3))
    check("log28_less_than_ten_thirds", low > 28)
    check("flat_log_margin", 9-R(38, 15) > R(5, 4))
    pi_hi = R(22, 7)
    check("Fourier_graph_high", 120*pi_hi < 400)
    check("Fourier_graph_low", 4*R(5, 4) == 5)
    check("fourth_to_flat", 30*pi_hi/R("13.9") < 7)
    check("scale_variation_budget", 2*(R("1.18")-1) == R(".36"))
    check("principal_Weyl_budget", 4*R("1.12")/R("13.9") < R(".33"))
    check("joint_principal_budget", R(".36")+R(".33") < R(".7"))
    check("interior_chart_location", R(1, 4)-R("1.1")/128 > R(1, 8))
    check("probe_support_margin", R(1, 128)+R(1, 1024) < R(1, 64))
    # Retain the negative terms in the degree-six seed on 0<=u<=1/4.
    check("actual_seed_outer_polynomial_lower",
          1-2*R("2.002")/4**3-R("3e-15") > R(".9"))
    check("canonical_distance_integrand", R("6.1")/R(".9") < 3**2)
    check("radius_anchor_before_inner_region", R(3, 28) < R(1, 8))
    check("actual_Weyl_principal_interval", R(10, 3)+R(1, 6) < 4)
    # Counting estimate discards only the negative part of the derivative weight.
    E, t = S.symbols("E t", positive=True)
    zero("tower_second_exponential_moment", S.integrate(4*t*E**4*S.exp(-t*E), (E, 0, S.oo))
         -96/t**4)
    check("Bessel_derivative_tower_reserve", R(2, 9)*96/4 < 1000)
    # Quantitative coordinate kernel comparison, including both density factors.
    distortion = (2*R("1.1")*100*1000/R(".9")**7
                  +2*100*72/R(".9")**6
                  +4*R("1.1")*72/R(".9")**4)
    check("fixed_chart_gain_one_kernel", distortion < 10**12, distortion)


def probe_reserves():
    x, y, t = S.symbols("x y t")
    indices = [(i, n-i) for n in range(4) for i in range(n+1)]
    probes = [x**i*y**j/(S.factorial(i)*S.factorial(j)) for i, j in indices]
    matrix = S.Matrix([[S.diff(p, x, i, y, j).subs({x: 0, y: 0})
                        for p in probes] for i, j in indices])
    check("ten_coordinate_jet_probes", matrix == S.eye(10))
    d = [11, 220, 3960, 63360, 887040, 10644480]
    radial = [1, 2, 6, 36, 300, 4000]
    cut = [S.Integer(1)]
    for j in range(1, 7):
        cut.append(2**j*sum(d[k-1]*S.bell(j, k, tuple(radial[:j-k+1]))
                            for k in range(1, j+1)))
    for j in range(7):
        # A degree<=3 normalized monomial has every derivative coefficient <=3!.
        bound = 6*sum(S.binomial(j, ell)*1024**(j-ell)*cut[j-ell]
                      for ell in range(min(j, 3)+1))
        check("rescaled_probe_jet_"+str(j), bound < 10**32, bound)
    yy = 2*t
    for degree in range(5):
        rhs = 400*(1-128*yy)**-2*S.diff(yy, t)**2
        coeff = S.expand(S.series(rhs, t, 0, degree+1).removeO()).coeff(t, degree)
        yy += coeff*t**(degree+2)/((degree+1)*(degree+2))
    jets = [S.diff(yy, t, j).subs(t, 0) for j in range(1, 7)]
    for n in [4, 6]:
        bound = 10**32*sum(2**k*S.bell(n, k, tuple(jets[:n-k+1])) for k in range(1, n+1))
        check("complete_flow_probe_X"+str(n), bound < 10**56, bound)
    check("model_probe_point_action", (1+R("2e-5"))*10**56 < 2*10**56)
    check("frozen_probe_point_action", (1+R("2e-5"))*2**12*10**32 < 2*10**56)
    # Covariant derivatives through three, metric density and the four seminorms.
    gamma = 100
    dgamma = 100*2*128
    conversion = 16*(1+1+(1+4*gamma)+(1+16*gamma+16*dgamma+64*gamma**2))
    check("physical_H3_chart_conversion", conversion < 10**12)


def strong_aggregation():
    pi_hi = R(22, 7)
    r0 = R("1e-20")
    near = 2*pi_hi*10**12*r0/24
    check("near_fourth_graph_absorption", 7*near < R("1e-6"))
    check("near_probe_point_reserve", near*4*10**32 < 10**26)
    tail_geo = 2*pi_hi*18/(R(".9")*r0)**4
    tail_flat = 2*pi_hi*1000/(4*r0**4)
    check("full_bounded_tail_Schur", 4*(tail_geo+tail_flat) < 10**86)
    check("all_ten_contact_coefficients", 4*10**56+10**26+10**86 < 10**87)
    check("model_lower_graph_budget", 10*10**87+10**86 < 10**90)
    check("actual_lower_graph_budget", 10**90+2*10**48*10**12 < 10**91)
    check("flat_outside_observation_tail", 36*pi_hi*256**4 < 10**13)
    check("principal_absorption_remaining", R(".36")+R(".33")+R("1e-6") < R(3, 4))
    check("flat_graph_lower_budget", 4*(10**91+10**13) < 10**93)
    check("fourth_log_intermediate_lower", 400*10**93+5 < 10**96)
    check("fourth_log_intermediate_forcing", 400*4 == 1600)
    eps = R(1, 2*10**96)
    check("interpolation_half_absorption", 10**96*eps == R(1, 2))
    check("final_value_reserve", 2*10**96*(2/eps)**3 < 10**390)
    check("final_forcing_reserve", 2*1600 < 4000)
    # Fourier interpolation in the two regions, encoded as exact monomial identities.
    q, ee, ell = S.symbols("q ee ell", positive=True)
    zero("high_frequency_interpolation_ratio", q**3/(ee*q**4*ell)-1/(ee*ell*q))
    zero("low_frequency_interpolation_threshold", q**3*((ee*ell*q)**-3)-(ee*ell)**-3)
    # Independent Taylor integral check through degree three.
    x, r, s = S.symbols("x r s")
    for n in range(13):
        p = x**n
        residual = p.subs(x, x+r)-sum(r**j*S.diff(p, x, j)/S.factorial(j) for j in range(4))
        integral = r**4*S.integrate((1-s)**3*S.diff(p, x, 4).subs(x, x+s*r), (s, 0, 1))/6
        zero("degree_three_Taylor_integral_"+str(n), residual-integral)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 compact massive logarithmic graph",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        flat_multiplier()
        continuum_and_geometry()
        probe_reserves()
        strong_aggregation()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 compact massive logarithmic graph",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "compact_massive_graph": good, "S166_complete": False, "goal_complete": False,
        "full_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["complete local TT operator and patch commutators",
                       "outer fast fourth/log estimate", "higher York reconstruction",
                       "full inverse/residual composition", "nonlinear signed compatibility",
                       "physical requirements"],
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
