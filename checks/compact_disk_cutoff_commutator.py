"""S166 one-patch cutoff bound: exact finite identities and numerical reserves."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "compact-disk-cutoff-commutator"
OWN = ["checks/compact_disk_cutoff_commutator.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("complete-compact-disk-graph", 217),
    ("strong-interior-quantum-remainder", 141),
]
PROTOCOL = {
    "scope": "one-patch localization of the complete fixed-area disk Hessian",
    "cutoff": "beta(6,6), inner 1/1024, outer 1/512",
    "larger_patch": "inner 1/128, outer 1/64",
    "near_radius": "1e-30",
    "far_input": "whole physical half L2(gbar), no far derivatives",
    "density_and_finite_contacts": "retained",
    "global_cover_absorption": False,
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
    expression = S.expand(expression)
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


def cutoffs_and_probes():
    t = S.symbols("t")
    poly = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    zero("cutoff_derivative_positive", S.diff(poly, t)-2772*t**5*(1-t)**5)
    zero("cutoff_zero", poly.subs(t, 0))
    zero("cutoff_one", poly.subs(t, 1)-1)
    bounds = [11, 220, 3960, 63360, 887040, 10644480]
    radial = [1, 2, 6, 36, 300, 4000]
    constants = [S.Integer(1)]
    for j, bound in enumerate(bounds, 1):
        pp = S.Poly(S.diff(poly, t, j), t)
        d = pp.degree()
        coeff = {i[0]: c for i, c in pp.terms()}
        bernstein = [sum(coeff.get(i, 0)*S.binomial(k, i)/S.binomial(d, i)
                         for i in range(k+1)) for k in range(d+1)]
        check("cutoff_Bernstein_"+str(j), max(abs(b) for b in bernstein) <= bound)
        if j < 6:
            zero("cutoff_left_join_"+str(j), pp.as_expr().subs(t, 0))
            zero("cutoff_right_join_"+str(j), pp.as_expr().subs(t, 1))
        constants.append(2**j*sum(bounds[k-1]*S.bell(j, k, tuple(radial[:j-k+1]))
                                 for k in range(1, j+1)))
    for j in range(7):
        check("chi_coordinate_jet_"+str(j), 512**j*constants[j] < 10**32)
        for degree in range(4):
            # Every normalized monomial component derivative is <=1 on its support.
            upper = sum(S.binomial(j, k)*1024**(j-k)*constants[j-k]
                        for k in range(min(j, degree)+1))
            check(f"probe_coordinate_jet_{degree}_{j}", upper < 10**32, upper)
    check("chi_Lipschitz", 2*22*512 < 10**6)
    check("chi_supported_in_V", R(1, 512) < R(1, 256))
    check("probe_support_in_fixed_chart", R(1, 128)+R(1, 1024) < R(1, 64))
    check("near_active_inside_eta", R(1, 512)+2*R("1e-30") < R(1, 128))
    x, y = S.symbols("x y")
    indices = [(i, j-i) for j in range(4) for i in range(j+1)]
    check("ten_jet_indices", len(indices) == 10)
    matrix = S.Matrix([[S.diff(x**a*y**b/(S.factorial(a)*S.factorial(b)), x, c, y, d)
                        .subs({x: 0, y: 0}) for a, b in indices] for c, d in indices])
    check("complete_two_dimensional_probe_matrix", matrix == S.eye(10))
    check("probe_product_old_chart", 2*2**5*10**64 < 10**70)
    for j in range(1, 6):
        check("chart_transition_jet_"+str(j),
              2*S.factorial(j-1)*512**(j-1) < S.factorial(j)*10000**j)
        jets = [S.factorial(i)*10000**i for i in range(1, j+1)]
        bell = sum(2**k*S.bell(j, k, tuple(jets[:j-k+1])) for k in range(1, j+1))
        zero("chain_Bell_identity_"+str(j),
             bell-S.factorial(j)*10000**j
             *sum(2**k*S.binomial(j-1, k-1)/S.factorial(k) for k in range(1, j+1)))
        check("chain_Bell_envelope_"+str(j), bell <= 2*3**(j-1)*S.factorial(j)*10000**j)
        check("recentered_probe_U_"+str(j), 10**70*bell < 10**100)
    # Physical distance of the probe support is below the recentered cap's core.
    check("probe_in_recentered_cap", R("1.1")**2/1024 < R(1, 128))
    check("transition_margin_old_chart",
          R(1, 128)+R(1, 1024)+R("1.1")**2/512 < R(1, 32))
    check("transition_margin_new_chart",
          R("1.1")**2/1024+R(1, 512) < R(1, 32))


def distribution_and_differential_identities():
    t = S.symbols("t")
    aa = S.symbols("a0:5")
    vv = S.symbols("v0:6")
    # The multiplier difference has a0=0. Its product's degree-four jet
    # depends on v0..v3 only, including for a general fourth contact.
    chi_difference = sum(aa[j]*t**j for j in range(1, 5))
    v = sum(vv[j]*t**j for j in range(6))
    fourth = S.expand(chi_difference*v).coeff(t, 4)
    zero("fourth_input_jet_cancels", S.diff(fourth, vv[4]))
    zero("fourth_product_all_lower_jets",
         fourth-(aa[1]*vv[3]+aa[2]*vv[2]+aa[3]*vv[1]+aa[4]*vv[0]))
    c, h, density, chi_value, q_value = S.symbols("c h density chi_value q_value")
    zero("density_contact_commutator", 2*chi_value*q_value*density
         -chi_value*2*q_value*density)
    zero("probe_direction_zero_value", (chi_value-chi_value)*q_value)
    x, y = S.symbols("x y")
    ch, field = S.Function("chi")(x, y), S.Function("v")(x, y)
    for order in range(1, 5):
        for nx in range(order+1):
            ny = order-nx
            actual = S.diff(ch*field, x, nx, y, ny)-ch*S.diff(field, x, nx, y, ny)
            expected = sum(S.binomial(nx, bx)*S.binomial(ny, by)
                           *S.diff(ch, x, bx, y, by)*S.diff(field, x, nx-bx, y, ny-by)
                           for bx in range(nx+1) for by in range(ny+1) if bx+by)
            zero(f"full_differential_commutator_{nx}_{ny}", actual-expected)
            check(f"commutator_input_order_{nx}_{ny}",
                  max(nx-bx+ny-by for bx in range(nx+1)
                      for by in range(ny+1) if bx+by) == order-1)
            check(f"commutator_binomial_sum_{nx}_{ny}",
                  sum(S.binomial(nx, bx)*S.binomial(ny, by)
                      for bx in range(nx+1) for by in range(ny+1) if bx+by)
                  == 2**order-1)
    r = S.symbols("r")
    for degree in range(11):
        f = x**degree
        remainder = f.subs(x, x+r)-sum(r**j*S.diff(f, x, j)/S.factorial(j) for j in range(4))
        integral = r**4*S.integrate((1-t)**3*S.diff(f, x, 4).subs(x, x+t*r), (t, 0, 1))/6
        zero("Taylor_integral_"+str(degree), remainder-integral)


def kernel_and_contact_budgets():
    pi_hi, rr = R(22, 7), R("1e-30")
    check("local_kernel_coordinate_conversion",
          4*(72/R(".9")**6+10**8*R(1, 1024)**2/R(".9")**4) < 10**12)
    check("local_commutator_gain", 10**12*10**6 == 10**18)
    tail = 100*(3*10**24+4*pi_hi*(72/(R(".9")*rr)**4+10**8/(R(".9")*rr)**2))
    check("full_tail_two_rows", tail < 10**30*rr**-4, tail)
    check("full_tail_numeric", 10**30*rr**-4 == 10**150)
    check("near_strong_with_point_reserve", 4*2*pi_hi*10**18/24 < 10**20)
    check("near_strong_numeric", 10**20*rr == R("1e-10"))
    check("actual_zero_value_probe_Cauchy", 6*10**22*10**100 == 6*10**122)
    check("probe_near_action", 10**20*rr*10**32 == 10**22)
    check("each_contact_coefficient", 6*10**122+10**22+10**150 < 10**152)
    check("ten_contacts_H3_reserve", 10*10**152 < 10**154)
    check("all_lower_local_commutators", 7*10**32*10**60 == 7*10**92)
    check("fourth_diagonal_coefficient_sum", 4*(R("1.12")/60+4*R("1.12")/90) < 1)
    check("fourth_local_commutator", 15*10**32 < 10**35)
    check("complete_local_commutator", 7*10**92+10**35 < 10**94)
    check("full_pair_commutator_H3", 10**154+10**94 < 10**155)
    xi, yi = S.symbols("xi yi")
    polynomial = S.Poly((1+xi**2+yi**2)**3, xi, yi)
    check("H3_derivative_coefficients", all(c > 0 for c in polynomial.coeffs())
          and sum(polynomial.coeffs()) == 27)
    check("H3_cutoff_multiplication", 27*(8*10**32)**2 < 10**80)


def graph_and_physical_budgets():
    check("graph_commutator_derivative", 2000*R("1e-10") == R("2e-7"))
    check("graph_commutator_H3", 2000*10**155+2*10**96*10**40 < 10**160)
    check("graph_commutator_value", 2000*10**150 == 2*10**153)
    eps = R("1e-167")
    check("interpolation_derivative", R("2e-7")+10**160*eps < R("1e-6"))
    check("interpolation_value", 2*10**160*(2/eps)**3+2*10**153 < 10**670)
    check("input_measure_conversion", R(5, 3)**2 < 4*3)
    check("physical_source_final", 2000*10**4 == 2*10**7)
    check("larger_patch_not_absorbed", PROTOCOL["global_cover_absorption"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 compact disk cutoff commutator",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        cutoffs_and_probes()
        distribution_and_differential_identities()
        kernel_and_contact_budgets()
        graph_and_physical_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 compact disk cutoff commutator",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "one_patch_cutoff_bound": good, "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["global compact cover absorption", "higher York maps",
                       "outer fast fourth/log estimate", "full inverse/residual composition",
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
