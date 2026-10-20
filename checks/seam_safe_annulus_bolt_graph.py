"""S166 seam-safe annulus/bolt graph: frozen arithmetic and full local Euler audit."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "seam-safe-annulus-bolt-graph"
OWN = ["checks/seam_safe_annulus_bolt_graph.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("outer-mass-column-source-bound", 244),
    ("compact-cover-graph-absorption", 366),
    ("strong-interior-quantum-remainder", 141),
    ("compact-massive-logarithmic-graph", 236),
    ("complete-compact-disk-graph", 217),
]
R = S.Rational
TINY = R(1, 10**170)
EPS = R(1, 10**160)
RAD = EPS/32
ROWS = []
PROTOCOL = {
    "scope": "canonical patch graph on Y>=r/4, including all intervening annulus and bolt",
    "fixed_scales": {"epsilon": "1e-160", "r": "epsilon/32", "chart_scale": "1e-170"},
    "interface_support": "r/8<Y<r/4, wholly inside chi_outer=1",
    "graph": "P <=1e8 F/C +1e-18000 B +1e381000 N",
    "actual_domain": "original paired resolvent; small caps only auxiliary",
    "Weyl_allowance": "abs(b)<400; input split at exp(2000)",
    "full_inverse": False,
    "nonlinear_contraction": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
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
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    found.update(["notes/general-interior-dirac-response.md",
                  "notes/curved-remote-dirac-response.md",
                  "checks/complete_compact_disk_graph.py"])
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


def charts_and_support():
    t, r = TINY, RAD
    check("one_sided_complex_domain", R(1, 16)+R(".02") < R(".1"))
    check("coordinate_derivative_near_identity", R(".002")/R(".999") < R(".003"))
    check("inverse_boundary_margin", R(".02")*(1-R(".003")) > R(".019"))
    check("log_ratio_bound", 2*R(".001")/(1-R(".001")) < R(".003"))
    check("normalizing_first_jet", R(".003")/R(".01") < 1)
    check("quadratic_chart_argument", t/32+(t/32)**2/2 < R(".01"))
    sig = R(".003")+t/32/(1-t/32)
    check("complex_normal_chart_exponential", sig/(1-sig) < R(".026"))
    check("interior_chart_before_seam", R("1.1")/32 < R(1, 16))
    check("small_chart_before_seam", R("1.1")*t/32 < r/8)
    check("eta_interface_lower", r/4-R("1.1")*t/64 > r/8)
    check("collar_interface_upper", r/4 < r)
    check("transition_fully_in_cover", r/4 <= 16*r)
    check("u_lower_beats_chart_scale", r/24 > t)
    # e^400 exceeds 24/r times e^(1/6); e^(1/6)<6/5.
    exp400_lower = sum(R(400)**k/S.factorial(k) for k in range(601))
    check("Weyl_bound_400", exp400_lower > R(6, 5)*24/r)
    check("exp_one_sixth", 1/(1-R(1, 6)) == R(6, 5))
    check("packing_count", 8/(R(3)*R(1, 4096)**2) < 10**8)
    dd = R(1, 4096)
    distance = R("1.1")*(R(1, 64)+R(1, 512))
    check("two_sided_neighbor_count", 8*((distance+dd/2)/(dd/2))**2 < 10**6)
    check("packing_chart_margin", R("1.1")*(distance+dd/2) < R(1, 32))
    check("finite_cover_size", 10**44*t**-2 == 10**384)
    check("square_root_cover_reserve", (10**192)**2 == 10**384)
    z, w, bb = S.symbols("z w b")
    # Independent complex variables, retaining both log Jacobian halves.
    h = lambda z: z-bb*z*z/2
    lin = bb*(h(z)+h(w))/2+(S.log(1-bb*z)+S.log(1-bb*w))/2
    zero("normal_chart_value", lin.subs({z: 0, w: 0}))
    zero("normal_chart_first_z", S.diff(lin, z).subs({z: 0, w: 0}))
    zero("normal_chart_first_w", S.diff(lin, w).subs({z: 0, w: 0}))


def scaled_first_cauchy():
    from functools import lru_cache
    scale = S.symbols("t", positive=True)
    p = {1: R(".02"), 2: R(9), 3: R(20000), 4: R(10**6)}
    b = {1: R(2), 2: R(100), 3: R(20000), 4: R(10**8), 5: R(10**12)}

    @lru_cache(None)
    def coeff(n, degree, scaled):
        ss = scale if scaled else S.Integer(1)
        if n == 1:
            return b[degree+1]*ss**(-degree-1)
        if degree == n-1:
            return (2/ss)**n
        return b[degree-n+2]*ss**(-degree+n-2)*(2/ss)**(n-1) \
            +sum(p[a]*ss**(-a)*coeff(n-1, degree-a, scaled)
                 for a in range(1, degree-n+2))

    for n, expected in enumerate([10**12, 204180000, 82644, R("808.0808")], 1):
        zero("Taylor_recursion_base_"+str(n), coeff(n, 4, False)-expected)
        zero("Taylor_recursion_scale_"+str(n), coeff(n, 4, True)-scale**-5*expected)
    rm = R("1.02")
    budget = 20*rm**2*coeff(1, 4, False) \
        +sum(28*S.factorial(5)*10**5*5**(n-2)*rm**(n+1)*coeff(n, 4, False)
             for n in range(2, 5)) \
        +28*S.factorial(5)*10**5*2**5*5**3*rm**6/(1-5*rm*R(".06"))
    check("full_Born_remainder_budget", budget < 8*10**16)
    check("Neumann_denominator", 5*rm*R(".06") < 1)
    # Total differential degree of every monomial in Q_phi is four.
    for name, orders in {"B4": [4], "A_squared": [2, 2], "p_dot_T": [1, 3],
                         "p_squared_A": [1, 1, 2]}.items():
        check("curvature_contact_scaling_"+name, sum(orders) == 4)
    check("small_evaluation_disk", R(1, 2560) < R(1, 256))
    check("Q_over_M", 1+R(256, 18) < 16)
    check("Caccioppoli_budget", 4*40**2+8 < 81**2)
    check("second_cutoff_budget", 5+2*80*81+16000 < 30000)
    check("resolvent_evaluation_exponential", 1/(1-R(1, 18)) < R(4, 3))
    check("small_scale_row", 4*10**6*16 < 10**8)
    check("actual_cap_source", 2*(10**8)**2*256*R("1.02") < 10**20)
    cauchy = 200*(R(2, 3)*10**17+10**20*24*256**3+10**17)
    check("complete_first_Cauchy_constant", cauchy < 10**40)
    for power in [4, 5, 6]:
        check("first_Cauchy_scale_reserve_"+str(power), power <= 10)
    m, v, source, dt = S.symbols("m v source dt")
    zero("fixed_area_Hessian_density", (dt+2*v*source).subs(v, 0)-dt)
    delta, x = S.symbols("delta x", positive=True)
    zero("angular_exponential_integral", S.integrate(delta*x**3*S.exp(-delta*x),
                                                     (x, 0, S.oo))-6/delta**3)


def massive_and_commutator_budgets():
    t, r = TINY, RAD
    first = 10**40*t**-10 * 10**18*t**-6
    geo = 2*10**46*t**-12
    near = 2*10**6*t*10**32*t**-6
    tail = 10**30*t**-4
    check("scaled_probe_first_Cauchy", first < 10**2790)
    check("scaled_probe_geodesic", geo < 10**2100)
    check("five_actual_contacts", first+geo+near+tail < 10**2800)
    pzero = 10**20*400**2*(1+8/r)
    check("constant_contact_retained", pzero < 10**200)
    cov = pzero+5*10**2800+2*10**6*t+tail
    check("strong_remainder_coordinate_conversion", cov*2*10**12*t**-6 < 10**4000)
    rr = R(1, 10**1000)
    check("model_split_inside_chart", rr < t/1024)
    check("model_freezing_Taylor_gain", R(2*22, 7*24)*10**12/t*rr*7 < R(1, 10**6))
    check("model_full_tail", 10**30*rr**-4 == 10**4030)
    probe_model = 4*10**56*t**-12
    probe_near = 10**14/t*rr*10**32*t**-6
    check("all_ten_model_contacts", probe_model+probe_near+10**4030 < 10**4040)
    check("complete_model_lower", 10*10**4040+10**4030 < 10**4050)
    check("flat_outside_source", 10**13*t**-4 < 10**4050)
    # e^.108 < 1/(1-.108); no floating point exponential enters.
    exp108 = 1/(1-R(".108"))
    check("large_Weyl_high_frequency", 400*exp108/(2000+9-R(38, 15)) < R(".23"))
    # e^2 <= 1+2+2 + tail starting at k=3 with ratio <=1/2 <10.
    check("e_squared_below_ten", 1+2+2+R(8, 6)/(1-R(1, 2)) < 10)
    check("Weyl_low_frequency_prefactor", 400*exp108/(30*3) < 10**10)
    check("principal_absorption", R(".36")+R(".23")+R(1, 10**6) < R(3, 4))
    check("supported_full_graph_lower", 400*4*(10**4050+10**4000+10**4010)
          +1600*10**2000+200*10**2000+10 < 10**4120)
    comm_r = R(1, 10**20000)
    check("commutator_split_inside_chart", comm_r < t/2048)
    check("commutator_near_gain", 10**20/t*comm_r < R(1, 10**19800))
    check("commutator_full_tail", 10**30*comm_r**-4 == 10**80030)
    actual_point = 10**140*t**-22
    check("commutator_actual_point", actual_point < 10**4000)
    check("commutator_near_probe", 10**20/t*comm_r*10**32*t**-6 < 1)
    check("commutator_all_contacts", actual_point+1+10**80030 < 10**80040)
    local_comm = 16*10**2000*10**32*t**-6
    check("complete_commutator_H3", 10*10**80040+local_comm < 10**80050)
    check("cutoff_H3_multiplier", 10**32*t**-6 < 10**1100)
    check("complete_one_patch_H3", 2000*10**80050+10**4120*10**1100 < 10**80060)
    check("complete_one_patch_near", 2000*10**20/t*comm_r < R(1, 10**19700))
    h = R(1, 10**100000)
    check("one_patch_interpolated_derivative", R(1, 10**19700)+10**80060*h
          < R(1, 10**19000))
    check("one_patch_interpolated_value", 2*10**80060*(2/h)**3+10**80040 < 10**380100)
    # Explicit homogeneous derivative ledger for the commutator probe.
    for j in range(1, 6):
        product = 2**j*10**64
        transition = 2*3**(j-1)*S.factorial(j)*10**(4*j)
        check("commutator_probe_U_constant_"+str(j), product*transition*512 < 10**100)
        check("commutator_probe_U_scale_"+str(j), 2*j+1 <= 12)


def complete_local_euler():
    spec = importlib.util.spec_from_file_location("previous_disk_graph",
                                                 ROOT/"checks/complete_compact_disk_graph.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    class Jets(module.Jets):
        def norm(self, expression, max_order, allocation):
            poly = S.Poly(S.expand(expression), *sorted(expression.free_symbols, key=str))
            symbols = poly.gens
            weights = {self.E: 2, self.Z: TINY**-2, self.e0: 1, self.k0: 60000,
                       self.ell: 1, self.charge: 2}
            total, highest, highest_background, max_z = S.Integer(0), 0, 0, 0
            for powers, coefficient in poly.terms():
                value, degree, background, zp, op = abs(coefficient), 0, 0, 0, 0
                for symbol, exponent in zip(symbols, powers):
                    if not exponent:
                        continue
                    if symbol in self.reverse:
                        field, idx = self.reverse[symbol]
                        order = sum(idx)
                        if field in ["theta", "T"]:
                            degree += exponent
                            highest = max(highest, order)
                            weight = 1
                        else:
                            background += order*exponent
                            weight = S.factorial(order)*(10**6/TINY)**order if order else \
                                (400 if field == "omega" else 1)
                            if field == "omega" and not order:
                                op += exponent
                    else:
                        weight = weights[symbol]
                        if symbol == self.Z:
                            zp += exponent
                    value *= weight**exponent
                assert degree == 1
                assert op <= 1
                highest_background = max(highest_background, background)
                max_z = max(max_z, zp)
                total += value
            check(allocation+"_input_order", highest <= max_order, highest)
            check(allocation+"_background_order", highest_background <= (4 if max_order == 3 else 2),
                  highest_background)
            check(allocation+"_exponential_degree", max_z <= (0 if max_order == 3 else 2), max_z)
            return total

    j = Jets()
    E, Z, e0, k0, ell, charge = j.E, j.Z, j.e0, j.k0, j.ell, j.charge
    theta, tt, omega, sigma = [j.f(f) for f in ["theta", "T", "omega", "sigma"]]
    ss = theta-tt/3
    rr, uu, vv = -2*E*j.lap(sigma), E*j.lap(omega), j.dot(omega, omega)
    dd, ff, ww, gg = E*j.lap(tt), E*j.lap(theta), j.dot(omega, ss), j.dot(ss, ss)
    dcoef = ss*(-rr/15+e0/9)+tt*((6*omega+1)*rr/45+11*vv/135)-11*ww/135
    lcoef = tt*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
    jcoef = rr/30-7*(uu+vv)/90-e0/18
    wcoef = -2*rr/15+14*(uu+vv)/45
    zcoef = tt**2*(omega*rr**2/30+rr*(uu+vv)/15-uu**2/5
                    -7*uu*vv/45-7*vv**2/90+2*k0*omega)+ss*tt*(-rr**2/30+2*k0)
    covq = ff**2/20-(omega+R(1, 6))*dd**2/30 \
        -(dd*dcoef+ff*lcoef+jcoef*gg+wcoef*tt*ww-R(7, 45)*ww**2+zcoef)/2
    density = S.expand(covq/E)
    operators = [S.expand(E*j.euler(density, field)/2) for field in ["theta", "T"]]
    leads = [E**2*j.lap(j.lap(theta))/20, -E**2*(omega+R(1, 6))*j.lap(j.lap(tt))/30]
    total = S.Integer(0)
    for field, op, lead in zip(["theta", "T"], operators, leads):
        bound = j.norm(S.expand(op-lead), 3, "Weyl_"+field)
        check("Weyl_reweighted_sum_"+field, bound < 10**1000)
        total += bound
        for idx, multiplicity in [((4, 0), 1), ((2, 2), 2), ((0, 4), 1)]:
            coefficient = E**2/20 if field == "theta" else -E**2*(omega+R(1, 6))/30
            zero("fourth_diagonal_"+field+str(idx), S.diff(op, j.v[field, idx])-multiplicity*coefficient)
            other = "T" if field == "theta" else "theta"
            zero("fourth_cross_"+field+str(idx), S.diff(op, j.v[other, idx]))
    qc = Z*(4*ss**2*(rr+6*vv)+48*ss*j.dot(omega, ss)
             +16*ss*j.dot(omega, tt)+12*j.dot(ss, ss)+8*j.dot(ss, tt)) \
        +8*e0*Z*(ss+tt)**2+6*ell*Z**2*(4*ss+2*tt)**2-8*charge*tt**2
    for field in ["theta", "T"]:
        op = S.expand(-E*j.euler(S.expand(qc/E), field)/8)
        bound = j.norm(op, 2, "classical_"+field)*2*10**33
        check("classical_reweighted_sum_"+field, bound < 10**1500)
        total += bound
    check("complete_local_coefficient_sum", total+100 < 10**2000)
    uj = [4, 2*10**4/TINY, 10**8/TINY**2, 2*10**11/TINY**3]
    for order in range(1, 5):
        bound = sum(S.factorial(k-1)*TINY**-k*S.bell(order, k, tuple(uj[:order-k+1]))
                    for k in range(1, order+1))
        check("omega_scaled_jet_"+str(order), bound < S.factorial(order)*(10**6/TINY)**order)


def chart_transfer_and_absorption():
    t = TINY
    # Retain physical z coordinates; scaling is in derivative constants only.
    native, trans, composed = R(10**7)/t, R(10**4)/t, R(10**12)/t**2
    hh, ee, qq = R(2*10**18)/t**2, R(10**28)/t**2, R(10**6)
    check("composed_cutoff_scale", 3*native*trans < composed)
    check("partition_derivative_scale", 2*10**6*composed == hh)
    check("localized_coefficient_scale", 2*composed+2*qq*hh < ee)
    for order in range(1, 6):
        base = tuple(S.factorial(k)*10**(4*k) for k in range(1, order+1))
        for k in range(1, order+1):
            zero(f"Bell_transition_identity_{order}_{k}", S.bell(order, k, base[:order-k+1])
                 -S.factorial(order)*S.binomial(order-1, k-1)*10**(4*order)/S.factorial(k))
    bvals = [R(2)]
    for k in range(1, 6):
        bvals.append(2*sum(2**j*S.bell(k, j, tuple(S.factorial(n)*trans**n
                                     for n in range(1, k-j+2))) for j in range(1, k+1)))
    for sob in [3, 5]:
        bounds = []
        for j in range(sob+1):
            bounds.append(sum(S.binomial(j, k)*4*S.factorial(j-k)*ee**(j-k)*bvals[k]
                              for k in range(j+1)))
        check("localized_H"+str(sob)+"_budget", 3**sob*max(bounds) < 10**160*t**-10)
    n = S.symbols("n", integer=True, positive=True)
    check("dyadic_log_schur", 32*(1+2*S.summation((n+1)/2**n, (n, 1, S.oo))) == 224)
    check("log_norm_transfer", 32*4*224 == 28672 and 28672 < 10**6)
    check("localized_log_reserve", 10**6*10**160*t**-10 < 10**1870)
    check("adjacency_reserve", 10**6*10**1870 == 10**1876)
    gain = R(1, 10**19000)
    absorption = gain*10**1876
    check("cover_absorption", absorption < R(1, 10**17000))
    check("cover_forcing", 2000*2000/(1-absorption) < 10**8)
    check("cover_interface", gain/(1-absorption) < R(1, 10**18000))
    check("cover_value", 10**380100*10**192/(1-absorption) < 10**381000)
    check("physical_source_conversion", 10*(24/RAD)**2 < 10**330)
    check("physical_forcing_coefficient", 10**8*10**330 == 10**338)
    check("compatible_value_substitution", R("5.58e9") < 10**10)
    a, ss = S.symbols("a s", nonzero=True)
    zero("partition_exact_identity", 1-a+a/ss*ss-1)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 seam-safe annulus/bolt graph",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        charts_and_support()
        scaled_first_cauchy()
        massive_and_commutator_budgets()
        complete_local_euler()
        chart_transfer_and_absorption()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S166 seam-safe annulus/bolt graph",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "annulus_bolt_graph": good, "interface_inside_outer_collar": good,
        "quantitative_interface_norm_transfer": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "S166_complete": False, "goal_complete": False,
        "unresolved": ["quantitative collar/patch norm transfers", "global graph absorption",
                       "higher reconstruction", "full inverse/residual composition",
                       "nonlinear signed compatibility", "physical requirements"],
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
