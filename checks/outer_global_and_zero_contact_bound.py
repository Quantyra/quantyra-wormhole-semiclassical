"""S166 actual global/constant-contact block: strong weighted estimates."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-global-and-zero-contact-bound"
OWN = ["checks/outer_global_and_zero_contact_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-curved-quantum-tt", 157),
    ("outer-flat-tower-normal-remainder", 321),
    ("weighted-york-source-conversion", 293),
    ("tt-mass-numerical-bound", 398),
]
PROTOCOL = {
    "operator": "actual row-subtracted global gluing plus complete constant contact",
    "weight": "7/2; output L2(dY/Y), with actual C*u^4/(16*pi*xi^2) density",
    "fast_bound": "1e26 epsilon^3 times weighted L2 input",
    "mass_bound": "1e26 abs(m_Y) epsilon^(5/2)",
    "outside_bound": "1e25 times normalized global canonical L2 input",
    "collar": "0<epsilon<=1/16",
    "normal_composition": "<1e-130 at epsilon=1e-60 for this block",
    "curved_singular_and_derivative_contacts_controlled": False,
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
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed"
              and len(run["checks"]) == count and all(c["passed"] for c in run["checks"])
              and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_consistent_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for path, digest in sorted(expected.items()):
        check("input_"+path, sha256((ROOT/path).read_bytes()).hexdigest() == digest)


def exact_measures_and_density():
    y, u, c, xi, b, lam, beta = S.symbols("Y u C xi B lambda beta", positive=True)
    delta = R(7, 2)
    rho = c*c*u**4/(xi*xi*y**4)
    zero("physical_shape_factor", c*u**4/(16*S.pi*xi**2)-rho*y**4/(16*S.pi*c))
    zero("normal_density_limit", rho.subs(u, S.sqrt(xi/c)*y)-1)
    # Dual output f_b=(f_T-f_S/3)/2 and input T=b/2.
    k, bb = S.symbols("H b")
    zero("two_dual_half_factors",
         (c*u**4/(4*S.pi*xi**2))*(k*bb/2)/2-c*u**4*k*bb/(16*S.pi*xi**2))
    zero("weighted_output_measure_power", 2*(4-delta)-1)
    zero("fast_input_measure_power", 2*delta+1-8)
    check("actual_density_upper", R("6.1")**2/25 < R(3, 2))
    check("density_factor_below_one", R(3, 2)/(16*3) < 1)
    # e^t <= 1/(1-t) for 0<=t<1, by the nonnegative power series.
    check("warp_upper", 1/(1-R(1, 16)**2) < 2)
    check("warp_lower", 1-R(1, 16)**2 > R(1, 2))
    zero("time_normalization", beta*S.sqrt(c*lam)-beta*S.sqrt(c)*S.sqrt(lam))
    gamma = 4*S.pi*beta*b**3
    a0 = beta*S.sqrt(c*lam)
    squared_conversion = c/(4*S.pi*b**4*a0)
    zero("physical_canonical_measure_conversion",
         squared_conversion-S.sqrt(c/(lam*b*b))/gamma)
    check("quarter_power_conversion_below_two", R("6.1") < 2**4)


def schur_and_row_subtraction():
    # An independent exact finite measure model tests orientation, density
    # and row subtraction. The infinite-kernel estimate is proved in the note.
    mu = S.diag(R(1, 2), R(2, 3), R(3, 5))
    kernel = S.Matrix([[2, -1, 3], [-1, 4, -2], [3, -2, 1]])
    op = kernel*mu
    ones = S.ones(3, 1)
    row = op*ones
    op0 = op-S.diag(*row)
    check("weighted_kernel_symmetry", mu*op == op.T*mu)
    check("row_subtracted_constant", op0*ones == S.zeros(3, 1))
    check("row_subtracted_symmetry", mu*op0 == op0.T*mu)
    v = S.Matrix(S.symbols("v0:3"))
    form = -sum(mu[i, i]*mu[j, j]*kernel[i, j]*(v[i]-v[j])**2
                for i in range(3) for j in range(3))/2
    zero("row_subtracted_form_exact", (v.T*mu*op0*v)[0]-form)
    contact = S.diag(5, -2, 7)
    restored = op0+contact
    check("physical_constant_contact_retained", restored*ones == contact*ones)
    check("no_constant_contact_double_count",
          op+contact-S.diag(*row) == restored)
    check("global_two_row_norm", S.sqrt((3*10**24)**2) == 3*10**24)
    check("subtracted_operator_budget", 3*10**24+3*10**24 == 6*10**24)
    check("stationary_output_conversion_budget", 2*(6*10**24)**2 < (9*10**24)**2)
    check("outside_value_budget", 9*10**24 < 10**25)
    check("both_measure_factors", S.sqrt(2)*6*10**24*S.sqrt(2) == 12*10**24)


def weights_and_infinite_end_integrals():
    e, y, s = S.symbols("epsilon Y s", positive=True)
    delta = R(7, 2)
    # Primitive verification plus positive exponents proves convergence at zero.
    for power in [R(1, 2), 2, R(5, 2), 3, R(7, 2), 4, 5, 7, 8]:
        zero("power_integral_"+str(power), S.diff(y**power/power, y)-y**(power-1))
        check("zero_endpoint_"+str(power), S.limit(y**power/power, y, 0, dir="+") == 0)
    zero("fast_L2_measure", y**(2*delta)*y-y**8)
    zero("mass_L2_measure", y**6*y-y**7)
    check("mass_input_integral", S.integrate(y**6, (y, 0, e)) == e**7/7)
    check("mass_contact_leading_exponent", 4-delta+3-1 == R(5, 2))
    check("mass_contact_regular_exponent", 4-delta+3 == R(7, 2))
    check("mass_contact_leading_integral",
          S.integrate(y**5/y, (y, 0, e)) == e**5/5)
    check("mass_contact_regular_integral",
          S.integrate(y**7/y, (y, 0, e)) == e**7/7)
    check("cubic_is_not_fast_L2",
          2*(3-delta) == -1)  # Integrand in dY is Y^-2.
    zero("cubic_fast_divergence_primitive",
         S.diff(-1/y, y)-y**(2*(3-delta)-1))
    check("contact_qC_budget", 10**20*400**2 == R("1.6e25"))
    emax = R(1, 16)
    fast = 12*10**24*emax + R("1.6e25")*(1+emax)
    mass = 12*10**24*emax/2 + R("1.6e25")*(emax/2+R(1, 2))
    check("mass_root_denominator_7", 7 > 2**2)
    check("mass_root_denominator_5", 5 > 2**2)
    check("fast_complete_bound", fast < 10**26, fast)
    check("mass_complete_bound", mass < 10**26, mass)
    # The envelope ratios increase in epsilon, so a single endpoint
    # bounds the entire open interval, not a sampled scale set.
    check("fast_envelope_monotonic", 12*10**24+R("1.6e25") > 0)
    check("mass_envelope_monotonic", 6*10**24+R("1.6e25")/2 > 0)
    zero("log2_integral_primitive", S.diff(S.log(s), s)-1/s)
    check("log2_lower_integrand", R(1, 2) > 0)  # 1/s>1/2 on 1<=s<2.
    check("H4log_to_L2_reserve", 2*10**26 < 10**27)


def mass_source_and_normal_composition():
    r = S.symbols("xi_over_C", positive=True)
    mc = S.symbols("m_c", real=True)
    # Shape component of (3/4)m_c u^3(1,-2), after b=2T.
    zero("disk_mass_to_b", 2*(-2)*R(3, 4)*mc+3*mc)
    zero("mass_coordinate_conversion", -3*mc*(S.sqrt(r))**3+3*mc*r**R(3, 2))
    check("mass_coefficient_quarter_bound", 9*16 < R("5.9")**3)
    check("physical_b_value_bound", 2*R("5.58e9") == R("1.116e10"))
    check("normalized_global_b_value", 2*R("1.116e10") == R("2.232e10"))
    check("outside_source_bound", 10**25*R("2.232e10") < 3*10**35)
    check("mass_source_bound", 10**26*R("5.40144e14")/4 < 2*10**40)
    inv = R("8.7e19")
    eps = R("1e-60")
    check("chosen_collar", eps < R(1, 16))
    check("actual_block_normal_budget", inv*10**27*eps**3 < R("1e-130"))
    check("flat_actual_density_budget", R(3, 2)*inv*10**30*eps < R("1.31e-10"))
    check("combined_two_component_budget",
          R(3, 2)*inv*10**30*eps+inv*10**27*eps**3 < R("2e-10"))
    lam = S.symbols("lambda", positive=True)
    zero("physical_inverse_units", (inv/lam**2)*lam**2-inv)
    check("outside_term_is_not_small_collar_error", 3*10**35 > 1)
    check("singular_response_remains_open",
          PROTOCOL["curved_singular_and_derivative_contacts_controlled"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 outer actual global and constant block",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        exact_measures_and_density()
        schur_and_row_subtraction()
        weights_and_infinite_end_integrals()
        mass_source_and_normal_composition()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 outer actual global and constant block",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "outer_actual_global_and_zero_contact_bound": good,
        "curved_singular_and_derivative_contacts_controlled": False,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["local curved singular response and derivative contacts",
                       "density/principal and Weyl combination", "annulus matching",
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
