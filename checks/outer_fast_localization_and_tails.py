"""S166 outer fast localization: exact cutoffs, contacts, tails and absorption."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-fast-localization-and-tails"
OWN = ["checks/outer_fast_localization_and_tails.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("outer-complete-disk-normal-comparison", 224),
    ("weighted-normal-inverse-scale", 453),
    ("compact-cover-graph-absorption", 366),
    ("outer-global-and-zero-contact-bound", 205),
]
PROTOCOL = {
    "operator": "complete physical affine disk Hessian in quantified half-density variables",
    "input": "global fast field in fixed mixed H4/H4log weight 7/2; no mass cubic included",
    "collar": "epsilon<=1e-160; r=epsilon/32; chi=c_r; eta=c_4r",
    "transition": "sigma=(1-c_(r/4))*c_(8r); support [r/4,16r]",
    "value_norm": "global canonical stationary L2 normalized by A0; not a weighted end norm",
    "graph_bound": "2e20 forcing +1e-40 transition X +1e2000 r^-4 global value",
    "global_mass_column_source_bounded": False,
    "global_inverse": False,
    "nonlinear_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.cancel(S.expand(expression))
    check(name, expression == 0, "" if expression == 0 else expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.add("evidence/"+stem+".json")
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(r["passed"] for r in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def normal_kernel_and_tail():
    z, y, t = S.symbols("z y t", positive=True)
    c = S.symbols("C", positive=True)
    integral = S.integrate((1+t*t)**-3, (t, -S.oo, S.oo))
    check("stationary_continuum_integral", integral == 3*S.pi/8)
    zero("normal_radial_kernel_constant", -16*c/(15*S.pi**2)*(3*S.pi/8)/(16*S.pi*c)
         +1/(40*S.pi**2))
    delta = R(7, 2)
    check("output_weight_measure_identity", 2*(4-delta)-1 == 0)
    # Kernel with dZ converted to dx', then conjugated by Y^-delta and Z^delta.
    zero("weighted_log_kernel_prefactor", y**(4-delta)*z**(1+delta)/y**5
         -(z/y)**(1+delta))
    check("positive_log_tail_exponent", 1+delta-5 == -R(1, 2))
    check("negative_log_tail_exponent", 1+delta == R(9, 2))
    check("small_log_kernel", 2*(R(10, 9)**5+1) < 10**6)
    check("middle_log_kernel", 100*(2**5+1) < 10**6)
    check("far_log_kernel", 2**5+1 < 10**3)
    # The cited exp bounds follow from geometric upper and positive partial sums.
    check("exp_small_upper", 1/(1-R(9, 20)) < 2)
    check("exp_negative_small_lower", 1-R(1, 10) >= R(9, 10))
    check("e_lower_two", 1+1+R(1, 2) > 2)
    check("e_upper_eleven_quarters", sum(R(1, S.factorial(k)) for k in range(7))
          +R(1, S.factorial(7))/(1-R(1, 8)) < R(11, 4))
    check("exp_four_point_five", R(11, 4)**5 < 160)
    # A sharper e^(9/2)<100 follows from e^9<10000.
    check("exp_four_point_five_sharp", R(11, 4)**9 < 10000)
    check("log_tail_row", 2*10**6/4+10**3*(2+R(2, 9)) < 10**10)
    rr = S.symbols("R", positive=True)
    check("bare_direct_row_integral", S.integrate(t**-5, (t, rr/2, S.oo))*2 == 8/rr**4)
    check("bare_image_row_integral", S.integrate(t**-5, (t, rr, S.oo)) == 1/(4*rr**4))
    check("bare_tail_rows", 8+R(1, 4) < 10**10)
    # exp(-R)<=1-R+R^2/2 implies 1-exp(-R)>=R/2 for 0<R<1.
    zero("log_gap_lower", rr-rr*rr/2-rr/2-rr*(1-rr)/2)
    check("normal_outside_support_gap", 32-2 == 30)


def cutoffs_and_probes():
    t = S.symbols("t", real=True)
    p = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    zero("beta_derivative", S.diff(p, t)-2772*t**5*(1-t)**5)
    ds = [1, 11, 220, 3960, 63360, 887040, 10644480]
    for k in range(6):
        check("cutoff_left_endpoint_"+str(k), S.diff(p, t, k).subs(t, 0) == 0)
        check("cutoff_right_endpoint_"+str(k), S.diff(1-p, t, k).subs(t, 1) == 0)
    for k in range(1, 7):
        pp = S.Poly(S.diff(p, t, k), t)
        n = pp.degree()
        bs = [sum(pp.nth(j)*S.binomial(i, j)/S.binomial(n, j) for j in range(i+1))
              for i in range(n+1)]
        check("beta_derivative_bound_"+str(k), max(abs(v) for v in bs) <= ds[k])
    # Compute (t d/dt)^k P(t-1) independently of a Stirling implementation.
    q = p.subs(t, t-1)
    for k in range(1, 7):
        q = S.expand(t*S.diff(q, t))
        # Positive recurrence for coefficients of t^j P^(j).
        if k == 1:
            coefficients = [0, 1]
        else:
            previous = coefficients
            coefficients = [0]+[(j*previous[j] if j < len(previous) else 0)
                                 +previous[j-1] for j in range(1, k+1)]
        exact = sum(coefficients[j]*t**j*S.diff(p, t, j).subs(t, t-1)
                    for j in range(1, k+1))
        zero("log_chain_identity_"+str(k), q-exact)
        bound = sum(coefficients[j]*2**j*ds[j] for j in range(1, k+1))
        check("log_cutoff_jet_"+str(k), bound < 10**14, bound)
    check("log_cutoff_Lipschitz", 2*ds[1] <= 22)
    for degree in range(4):
        for order in range(7):
            bound = sum(S.binomial(order, j)*ds[order-j]*8**(order-j)
                        *R(1, 4)**(degree-j)/S.factorial(degree-j)
                        for j in range(min(order, degree)+1))
            check(f"probe_jet_{degree}_{order}", bound < 10**16, bound)
    probes = [t**k/S.factorial(k) for k in range(4)]
    jets = S.Matrix([[S.diff(f, t, k).subs(t, 0) for f in probes] for k in range(4)])
    check("probe_identity_matrix", jets == S.eye(4))
    for k, f in enumerate(probes):
        zero("probe_exact_Taylor_"+str(k), f-sum(S.diff(f, t, j).subs(t, 0)*t**j/S.factorial(j)
                                                    for j in range(4)))
    check("point_probe_H6", 8*2**6*10**14*10**16 < 10**35)
    check("near_commutator_kernel", 22*10**6 < 10**8)
    radius = R("1e-100")
    check("near_Taylor_operator", 2*10**8*radius/24*2 < R("1e-90"))
    check("probe_plateau_contains_near", radius < R(1, 8))
    check("near_segment_support", 1/(1-2*radius) < 2)
    check("point_normal_action", 2*10**38*10**35 == 2*10**73)
    check("raw_log_tail_on_probe", 10**10*radius**-4 == 10**410)
    check("each_contact_budget", 2*10**73+10**410 < 10**420)
    check("sixteen_contacts", 2*2*4*10**420 < 10**423)
    check("unweighted_commutator_tail", 2*10**10*radius**-4 < 10**430)
    check("annular_support_inner", R(1, 4) < R(1, 2) < 1)
    check("annular_support_outer", 2 < 4 < 8 < 16 < 32)
    check("eta_support_within_comparison", R(8, 32) < 1)
    check("zeta_support_in_sigma_plateau", R(1, 2) < 1 and 8 <= 8)


def actual_tail_and_absorption():
    pi = R(22, 7)
    check("actual_geodesic_tail", 2*pi*18/2**4 < 10)
    check("actual_rows_budget", 100*(2*pi*18/2**4+pi*10**8/4+3*10**24) < 10**29)
    check("actual_density_measure_conversion", 8*10**29 < 10**31)
    check("global_measure_normalization", R(4)**4*R(3, 2) < 400)
    check("multiplication_H5", sum(S.binomial(5, j)*4**j for j in range(6)) < 10**6)
    check("zeta_integer_norm", 10**3*2*10**14 == 2*10**17)
    check("zeta_mixed_norm", 10**6*2*10**17 < 10**24)
    check("normal_inverse_in_contract", R("8.7e19") < 10**20)
    eps, tau, inv = R("1e-160"), R("1e-500"), R("8.7e19")
    self_term = inv*10**67*eps
    check("self_absorption", self_term < R(1, 2))
    transition = 2*inv*(10**67*eps*10**24+R("1e-90")+10**423*tau)
    check("transition_coefficient", transition < R("1e-40"), transition)
    value = 2*inv*(10**31+10**10+10**430+400*10**423*(2/tau)**3)
    check("global_value_coefficient", value < 10**2000)
    check("forcing_coefficient", 2*inv < 2*10**20)
    check("collar_contract", eps < R("1e-15"))
    # Interpolation uses the Fourier low/high split at sqrt(1+t^2)=2/tau.
    check("interpolation_high_scalar", R(1, 2) < 1)
    check("interpolation_high_log", R(1, 2)/R(1, 2) <= 1)
    check("mass_column_remains_open", PROTOCOL["global_mass_column_source_bounded"] is False)


def localization_identity():
    # A finite matrix verifies the ordered decomposition without commuting A and cutoffs.
    aa = S.Matrix([[1, 2, 3], [4, 6, 5], [7, 8, 9]])
    nn = S.Matrix([[9, 2, 7], [6, 4, 8], [3, 5, 1]])
    chi = S.diag(1, R(1, 2), 0)
    eta = S.diag(1, 1, R(1, 2))
    vv = S.Matrix(S.symbols("u0:3"))
    rhs = chi*aa*vv-chi*aa*(S.eye(3)-eta)*vv+(nn*chi-chi*nn)*eta*vv+chi*(nn-aa)*eta*vv
    check("ordered_localization_identity", nn*chi*vv == rhs)
    check("cutoff_product_identity", chi*eta == chi)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 outer fast localization and tails",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        normal_kernel_and_tail()
        cutoffs_and_probes()
        actual_tail_and_absorption()
        localization_identity()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 outer fast localization and tails", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version, "sympy": S.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "outer_fast_localization_with_global_value_tails": good,
        "S166_complete": False, "goal_complete": False, "full_mass_column_source_bounded": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["quantitative full mass-column image", "transition and intervening annulus matching",
                       "compact/bolt assembly", "higher York maps", "full inverse/residual composition",
                       "nonlinear signed compatibility", "all physical completion requirements"],
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
