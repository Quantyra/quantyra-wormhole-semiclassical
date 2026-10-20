"""S166 compact cover, logarithmic chart bounds and interface graph absorption."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "compact-cover-graph-absorption"
OWN = ["checks/compact_cover_graph_absorption.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("compact-disk-cutoff-commutator", 332),
    ("retained-benchmark-parameter-recovery", 25),
    ("weighted-york-source-conversion", 293),
]
PROTOCOL = {
    "compact_region": "Y>=1/4 including regular bolt",
    "cover": "any fixed maximal 1/4096-separated intrinsic net",
    "norm": "l2 of sum H4(theta) and H4log(T) on chi patches",
    "outer_interface": "eta_i zeta_0 w, supported in 7/32<Y<1/4",
    "overlap": "1000000, both neighbor rows",
    "near_radius": "1e-220",
    "log_mapping": "explicit orthogonal dyadic blocks using H3 and H5",
    "compact_derivative_absorption": True,
    "outer_interface_controlled": False,
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


def cover_geometry_and_period():
    d, pi_hi = R(1, 4096), R(22, 7)
    check("net_ball_in_chi_core", R("1.1")*d < R(1, 1024))
    check("small_ball_area_ratio", R(".9")**2/R("1.1")**2 > R(1, 2))
    check("small_ball_inside_half", R(1, 4)-d/2 > R(1, 8))
    check("area_upper_chart_ratio", R("1.1")**4 < 4)
    check("cover_count_per_area", 2/(3*(d/2)**2) < 10**8)
    distance = R("1.1")*(R(1, 64)+R(1, 512))
    check("neighbor_ball_in_chart", distance+d/2 < R(".9")/32)
    check("neighbor_rows_and_reverse_rows", 8*((distance+d/2)/(d/2))**2 < 10**6)
    check("interface_lower_distance", R(1, 4)-R("1.1")/64 > R(7, 32))
    check("source_support_in_known_region", R(1, 4)-R("1.1")/512 > R(1, 8))
    check("coordinate_margin_eta", R(1, 64)+R("1.1")**2/512 < R(1, 32))
    check("coordinate_margin_chi", R(1, 512)+R("1.1")**2/512 < R(1, 32))
    beta, cc, horizon, r = S.symbols("beta C R_h r", positive=True)
    zero("exact_canonical_area", S.integrate(beta*cc/r**2, (r, horizon, S.oo))-beta*cc/horizon)
    k0 = 2000
    check("benchmark_C0_upper", 1+12*pi_hi*k0**6 < 49*k0**6)
    check("benchmark_xi0_upper", R(7, 6) < 2)
    check("benchmark_sqrt_z_over_u", R("1.01") < 4*R(".99")**2)
    check("canonical_orbit_area", 16*pi_hi*k0**10 < 10**36)
    check("cover_count", 10**8*10**36 == 10**44)
    check("cover_count_square_root", (10**22)**2 == 10**44)
    # Exact cancellation of the physical length scale and of C0.
    q, eta, gnew, c0, b0, z, u, xi, kk = S.symbols(
        "q eta G C0 b0 z u xi K", positive=True)
    nu = 2*q*eta*gnew*c0/b0**3
    area_envelope = (8*S.pi/nu)*(c0*z/u)/(b0*S.sqrt(z))
    zero("area_benchmark_cancellation",
         area_envelope-4*S.pi*b0**2*S.sqrt(z)/(q*eta*gnew*u))
    zero("area_benchmark_power",
         (4*S.pi*b0**2/(q*eta*gnew)).subs(
             {b0: S.sqrt(gnew*kk**10*xi), q: kk**4, eta: 1/kk})
         -4*S.pi*kk**7*xi)


def cutoff_and_partition_jets():
    t = S.symbols("t")
    poly = 462*t**6-1980*t**7+3465*t**8-3080*t**9+1386*t**10-252*t**11
    d = [1, 11, 220, 3960, 63360, 887040, 10644480]
    radial = [1, 2, 6, 36, 300, 4000]
    const = [S.Integer(1)]
    for j in range(1, 6):
        pp = S.Poly(S.diff(poly, t, j), t)
        deg = pp.degree()
        coeff = {i[0]: c for i, c in pp.terms()}
        bernstein = [sum(coeff.get(i, 0)*S.binomial(k, i)/S.binomial(deg, i)
                         for i in range(k+1)) for k in range(deg+1)]
        check("beta_Bernstein_"+str(j), max(abs(x) for x in bernstein) <= d[j])
        const.append(2**j*sum(d[k]*S.bell(j, k, tuple(radial[:j-k+1]))
                             for k in range(1, j+1)))
        check("native_cutoff_factorial_jet_"+str(j),
              512**j*const[j] < S.factorial(j)*10**(7*j))
        for endpoint in [0, 1]:
            zero(f"beta_join_{j}_{endpoint}", S.diff(poly, t, j).subs(t, endpoint))
        transition = [S.factorial(i)*10**(4*i) for i in range(1, j+1)]
        composed = sum(2**k*S.factorial(k)*10**(7*k)
                       *S.bell(j, k, tuple(transition[:j-k+1])) for k in range(1, j+1))
        exact = S.factorial(j)*10**(4*j)*(2*10**7)*(1+2*10**7)**(j-1)
        zero("composed_cutoff_Bell_"+str(j), composed-exact)
        check("composed_cutoff_D_"+str(j), composed < S.factorial(j)*10**(12*j))
        check("sum_squares_jet_"+str(j),
              10**6*(j+1)*10**(12*j) <= (2*10**18)**j)
    q = 10**6
    for k in range(6):
        transition_bound = S.factorial(k)*4**(k+1)*sum(R(d[j], S.factorial(j)) for j in range(k+1))
        check("partition_scalar_transition_"+str(k),
              transition_bound <= 4*S.factorial(k)*q**k, transition_bound)
        check("partition_scalar_constant_region_"+str(k),
              S.factorial(k)*2**(k+1) <= 4*S.factorial(k)*q**k)
    ss, aa = S.symbols("s a")
    zero("partition_identity", (1-aa)+(aa/ss)*ss-1)
    zero("partition_on_compact_core", (1-aa).subs(aa, 1))
    zero("partition_reconstruction_factor", aa/ss*ss-aa)
    h = 2*10**18
    for j in range(1, 6):
        jets = [S.factorial(i)*h**i for i in range(1, j+1)]
        total = sum(4*S.factorial(k)*q**k*S.bell(j, k, tuple(jets[:j-k+1]))
                    for k in range(1, j+1))
        zero("partition_composition_exact_"+str(j),
             total-4*S.factorial(j)*h**j*q*(1+q)**(j-1))
        check("partition_composition_envelope_"+str(j),
              total <= 4*S.factorial(j)*(2*q*h)**j)
        check("three_factor_coefficient_jet_"+str(j),
              4*S.factorial(j)*(2*10**12+2*q*h)**j < 4*S.factorial(j)*10**(28*j))


def coordinate_operator_and_log_norm():
    b, ee = 10**4, 10**28
    chain = [S.Integer(2)]
    for k in range(1, 6):
        jets = [S.factorial(i)*b**i for i in range(1, k+1)]
        chain.append(2*sum(2**l*S.bell(k, l, tuple(jets[:k-l+1])) for l in range(1, k+1)))
    budgets = [sum(S.binomial(j, k)*4*S.factorial(j-k)*ee**(j-k)*chain[k]
                   for k in range(j+1)) for j in range(6)]
    for order in [3, 5]:
        check("integer_operator_norm_"+str(order),
              3**order*max(budgets[:order+1])**2 < 10**320)
    j, kk, nn, dd = S.symbols("j k n d", integer=True, nonnegative=True)
    # Universal dyadic inequalities, not sampled frequencies.
    zero("dyadic_log_ratio_forward", (1+dd)*(1+kk)-(1+kk+dd)-dd*kk)
    check("dyadic_log_ratio_reverse", S.ask(S.Q.nonnegative(dd)) is True)
    zero("log_lower_large_blocks", ((3+nn)-1)/2-((3+nn)+1)/4-nn/4)
    check("log_lower_block_zero", R(1, 2) > R(1, 4))
    check("log_lower_block_one", 1 > R(2, 4))
    check("log_lower_block_two", 1 > R(3, 4))
    zero("log_upper_all_blocks", 2*(j+1)-(j+2)-j)
    # log(2)>1/2 follows from exp(1/2)<sum (1/2)^n=2;
    # log(3)>1 follows from n!>=2^(n-1), with strictness for n>=3.
    check("log_two_geometric_majorant", S.summation(R(1, 2)**nn, (nn, 0, S.oo)) == 2)
    check("log_three_factorial_majorant",
          2+S.summation(R(1, 2)**(nn-1), (nn, 2, S.oo)) == 3)
    check("dyadic_low_radius", 3*4**0-1 >= 0)
    schur_log = 1+2*S.summation((nn+1)*R(1, 2)**nn, (nn, 1, S.oo))
    schur_plain = 1+2*S.summation(R(1, 2)**nn, (nn, 1, S.oo))
    check("infinite_log_Schur_sum", schur_log == 7)
    check("infinite_plain_Schur_sum", schur_plain == 3)
    check("H4log_dyadic_constant", 128*32*schur_log == 28672)
    check("H4_dyadic_constant", 16*32*schur_plain == 1536)
    check("common_dyadic_mapping_reserve", max(28672, 1536) < 10**6)
    check("complete_chart_mapping_reserve", 10**6*10**160 < 10**170)


def tightened_patch_and_absorption():
    near, eps = R("1e-220"), R("1e-1105")
    check("tightened_near_coefficient", 10**20*near == R("1e-200"))
    check("tightened_tail_coefficient", 10**30*near**-4 == 10**910)
    check("tightened_each_contact", 6*10**122+10**20*near*10**32+10**910 < 10**912)
    check("tightened_complete_H3_commutator", 10*10**912+10**94 < 10**915)
    check("tightened_patch_high", 2000*10**20*near == R("2e-197"))
    check("tightened_patch_H3", 2000*10**915+2*10**96*10**40 < 10**920)
    check("tightened_patch_value", 2000*10**910 == 2*10**913)
    check("tightened_interpolated_high", R("2e-197")+10**920*eps < R("1e-180"))
    check("tightened_interpolated_value", 2*10**920*(2/eps)**3+2*10**913 < 10**4240)
    check("two_neighbor_rows", 10**170*10**6 == 10**176)
    check("compact_high_absorption", R("1e-180")*10**176 == R("1e-4"))
    factor = 1/(1-R("1e-4"))
    check("absorption_factor", factor < 2)
    check("source_overlap_norm", 2*S.sqrt(10**6) == 2000)
    check("global_forcing_reserve", factor*2000*2000 < 10**8)
    check("interface_coefficient", factor*R("1e-180") < R("1e-179"))
    check("global_value_reserve", factor*10**4240*10**22 < 10**4270)
    check("physical_forcing_conversion", 10**8*10**4 == 10**12)
    check("compatible_value_substitution", 10**4270*R("5.58e9") < 10**4280)
    check("outer_interface_remains", PROTOCOL["outer_interface_controlled"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 compact cover graph absorption",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        cover_geometry_and_period()
        cutoff_and_partition_jets()
        coordinate_operator_and_log_norm()
        tightened_patch_and_absorption()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 compact cover graph absorption",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "compact_cover_absorption": good, "outer_interface_controlled": False,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["outer annulus and weighted fast end", "higher York maps",
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
