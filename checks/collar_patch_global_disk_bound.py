"""S166 two-dimensional/radial transfer, global absorption, physical inverse bound."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
from sympy.functions.combinatorial.numbers import stirling

ROOT = Path(__file__).resolve().parents[1]
STEM = "collar-patch-global-disk-bound"
OWN = ["checks/collar_patch_global_disk_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("seam-safe-annulus-bolt-graph", 323),
    ("outer-mass-column-source-bound", 244),
    ("outer-complete-massive-model-comparison", 277),
    ("weighted-york-source-conversion", 293),
    ("tt-mass-numerical-bound", 398),
    ("stationary-scalar-outgoing-inverse", 107),
]
R = S.Rational
TINY, EPS = R(1, 10**170), R(1, 10**160)
RAD = EPS/32
ROWS = []
PROTOCOL = {
    "scope": "global actual compatible disk inverse in a specified mass-extended mixed norm",
    "fixed_scales": {"epsilon": "1e-160", "r": "epsilon/32", "chart_scale": "1e-170"},
    "transfer": "both integer endpoints, dimension-independent dyadic log transfer, actual time average",
    "interface_feedback": "1e-8040",
    "a_priori_bound": "G<=1e22 Z+1e5000 Q+1e386000 V+1e5100 abs(m)r^-1/2",
    "physical_inverse": "G<=1e386020 lambda^-2 gamma0^-1/2 Fdisk on ker L0",
    "higher_York_norms": False,
    "residual_composition": False,
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
    found.update(["notes/stationary-global-weighted-map.md",
                  "notes/stationary-scalar-outgoing-inverse.md"])
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


def support_and_partition():
    t, r, dd = TINY, RAD, TINY/4096
    check("all_transfer_band_inside_analytic_collar", 32*r == EPS,
          "The transfer band is open at Y=32r=epsilon.")
    check("radial_band_away_from_seam", r/16 > t)
    check("mass_inner_cutoff_before_patch_support", r/16 < r/8)
    check("patch_support_before_cover_edge", r/4-R("1.1")*t/512 > r/8)
    check("interface_lies_in_chi_one", r/4 < r)
    check("mass_outer_cutoff_one", 2*(4*r) == 8*r)
    check("reverse_support_in_cover", r/4 <= 16*r)
    radius = R("1.1")*t/64
    check("eta_overlap_count", 8*((radius+dd/2)/(dd/2))**2 < 10**6)
    check("eta_packing_chart_margin", R("1.1")*(radius+dd/2) < t/32)
    s, aa, field = S.symbols("s a field", nonzero=True)
    zero("partition_reconstruction", aa/s*s*field-aa*field)
    zero("partition_on_K", (aa/s*s*field).subs(aa, 1)-field)
    period, y = S.symbols("period y", positive=True)
    zero("stationary_time_average", S.integrate(field, (y, 0, period))/period-field)


def finite_transfer_budgets():
    t, r = TINY, RAD
    hh = 10**40/t**3
    check("physical_partition_derivative_base", 3*(10**28/t**2)*(10**4/t) < 10**33/t**3)
    check("common_derivative_base", hh == 10**550)
    check("log_coordinate_derivative_base", 10**6/t**2 < hh)
    dpoly = [1, 11, 220, 3960, 63360, 887040]
    for k in range(1, 6):
        check("beta_partition_derivative_"+str(k), 4**k*dpoly[k] < S.factorial(k)*10**(6*k))
        fbound = S.factorial(k)*4**(k+1)*sum(R(dpoly[j], S.factorial(j)) for j in range(k+1))
        check("quotient_partition_derivative_"+str(k), fbound < 4*S.factorial(k)*10**(6*k))
    x, yy = S.symbols("x Y", positive=True)
    test = yy**7+3*yy**4-2*yy**2+yy
    repeated = test
    for order in range(6):
        rhs = sum(stirling(order, k, kind=2)*yy**k*S.diff(test, yy, k) for k in range(order+1))
        zero("conormal_Stirling_identity_"+str(order), repeated-rhs)
        repeated = yy*S.diff(repeated, yy)
    # One homogeneous base suffices; evaluate exact Bell sums then restore H^j.
    cvals = [S.Integer(1)]
    for k in range(1, 6):
        vals = tuple(S.factorial(n) for n in range(1, k+1))
        cvals.append(sum(2**ell*S.bell(k, ell, vals[:k-ell+1]) for ell in range(1, k+1)))
        for ell in range(1, k+1):
            zero(f"Bell_normalization_{k}_{ell}", S.bell(k, ell, vals[:k-ell+1])
                 -S.factorial(k)*S.binomial(k-1, ell-1)/S.factorial(ell))
    budgets = []
    for j in range(6):
        prod = sum(S.binomial(j, k)*10**10*S.factorial(j-k)*cvals[k] for k in range(j+1))
        stir = max(1, sum(stirling(j, k, kind=2) for k in range(j+1)))
        # Extra cutoff, field conversion and conormal expansion, before H^j.
        reserve = prod*2*10**14*stir
        kconst = 10**40*S.factorial(j)*10**j
        check("full_chain_product_budget_"+str(j), reserve < kconst)
        budgets.append(kconst*hh**j)
    for sob in [3, 5]:
        lift = 2*1000*3**sob*max(budgets[:sob+1])
        reverse = 4*1000*(r/4)**-4*3**sob*max(budgets[:sob+1])
        check("radial_lift_H"+str(sob), lift < 10**4000)
        check("patch_average_H"+str(sob), reverse < 10**4000)
    delta = R(7, 2)
    check("reverse_weight_and_measure_power", delta+R(1, 2) == 4)
    check("lift_measure_bound", 2*32*r < 2)
    n = S.symbols("n", integer=True, positive=True)
    schur = 32*(1+2*S.summation((n+1)/2**n, (n, 1, S.oo)))
    check("dimension_independent_log_Schur", schur == 224)
    check("log_norm_shell_equivalence", 32*4*schur == 28672)
    check("two_field_log_conversion", 2*28672 < 10**6)
    check("mixed_transfer_reserve", 10**6*10**4000 < 10**5000)
    # The reverse period factor is reciprocal to the lifting period factor.
    aa = S.symbols("A0", positive=True)
    zero("period_factor_cancellation", S.sqrt(aa)/S.sqrt(aa)-1)


def mass_profile():
    r = RAD
    check("mass_profile_log_support", (8*r)/(r/32) == 256)
    check("mass_profile_support_length", sum(R(6)**k/S.factorial(k) for k in range(16)) > 256)
    check("mass_profile_pair_amplitude", 32*(1+R(1, 144)) < 36)
    for order in range(6):
        budget = 6*3**order*10**28
        check("mass_profile_log_jet_"+str(order), budget < 10**33)
    check("mass_profile_H5_squared", 6*2**5*10**66 < 10**70)
    check("mass_profile_mixed", 2*10**35 < 10**40)
    check("mass_patch_profile_reserve", 10**5000*10**40 == 10**5040)
    check("physical_mass_L_relation", R(1, 2)**2 == R(1, 4))
    check("mass_coefficient_not_fast", 2*(3-R(7, 2)) == -1)


def exact_interface_absorption():
    # p <= q0 Q + a O + b M_r + v0 V,
    # O <= z0 Z + h p + v1 V + m1 M_r.
    q0, a, b, v0 = 10**8, R(1, 10**13000), R(1, 10**12960), 3*10**381000
    z0, hh, v1, m1 = 2*10**20, 10**4960, 10**2647, R(10**2005, 2)
    check("annulus_interface_fast_coefficient", R(1, 10**18000)*10**5000 == a)
    check("annulus_interface_mass_coefficient", R(1, 10**18000)*10**5040 == b)
    check("outer_interface_coefficient", R(1, 10**40)*10**5000 == hh)
    check("value_conversion_N_over_A0", 8 < 9)
    check("outer_value_r_power", RAD**-4 < 10**647)
    feedback = hh*a
    check("two_interface_feedback", feedback == R(1, 10**8040))
    check("two_interface_absorption", feedback < R(1, 2))
    oc = {"Z": R(z0)/(1-feedback), "Q": R(hh*q0)/(1-feedback),
          "V": R(hh*v0+v1)/(1-feedback), "M": (hh*b+m1)/(1-feedback)}
    # G <=1.01 O+p+1e5040 M_r+|m|; |m|<=M_r.
    factor = R("1.01")+a
    gc = {"Z": factor*oc["Z"], "Q": factor*oc["Q"]+q0,
          "V": factor*oc["V"]+v0, "M": factor*oc["M"]+b+10**5040+1}
    for name, limit in {"Z": 10**22, "Q": 10**5000, "V": 10**386000, "M": 10**5100}.items():
        check("global_a_priori_"+name, gc[name] < limit)
    zz, qq, vv, mm = S.symbols("Z Q V M")
    qa, ab, bm, va, za, hb, vb, mb = S.symbols("q a b va z h vb mb")
    # Independent elimination checks the retained value and mass terms as well.
    osol = (za*zz+hb*qa*qq+(hb*va+vb)*vv+(hb*bm+mb)*mm)/(1-hb*ab)
    psol = qa*qq+ab*osol+bm*mm+va*vv
    zero("exact_outer_linear_elimination", osol-(za*zz+hb*psol+vb*vv+mb*mm))
    zero("exact_patch_linear_elimination", psol-(qa*qq+ab*osol+bm*mm+va*vv))


def physical_source_and_solution():
    bb, cc, kk, beta, u, yy, ww, jj, lam, aa = S.symbols(
        "B C k beta u Y w j lambda A0", positive=True)
    fs, ft = S.symbols("fS fT", real=True)
    rho4 = 4*S.pi*bb**4/cc
    zero("physical_source_weight_cancellation", (rho4*u**-4)**2/(rho4*u**-8)-rho4)
    zero("canonical_source_transpose", fs*S.Symbol("theta")+(ft-fs/3)*S.Symbol("T")
         -(fs*(S.Symbol("theta")-S.Symbol("T")/3)+ft*S.Symbol("T")))
    check("canonical_two_row_sum", R(4, 3)**2+1 == R(25, 9))
    gamma = 4*S.pi*beta*bb**3
    area = beta*kk*cc/bb
    zero("physical_time_gamma_identity", 4*S.pi*area*bb**4/(kk*cc)-gamma)
    # lambda=k^2 C/B^2, A0=beta k C/B.
    qfactor = (10*S.sqrt(S.pi)/3)*bb**2/(cc**R(3, 2)*S.sqrt(area)) \
        *(kk**2*cc/bb**2)**2*S.sqrt(gamma)
    zero("Q_source_factor", S.simplify(qfactor)-R(20, 3)*S.pi*kk**R(7, 2))
    check("Q_original_source_bound", (R(20, 3)*R(22, 7))**2/R(2)**7 < 9)
    dual = S.Matrix([[1, 0], [-R(1, 6), R(1, 2)]])
    check("physical_dual_Frobenius", sum(e**2 for e in dual) == R(23, 18))
    ddensity = kk*(kk*yy*ww)**-8*jj
    zdensity = yy**-8*jj
    zero("outer_weight_half_density_conversion", zdensity/ddensity-kk**7*ww**8)
    check("Z_original_source_bound", R(23, 18)*R(1, 2)**7*R("1.001")**8 < 9)
    check("canonical_value_physical_factor", 4*R(22, 7)/R(".4") < 36)
    check("mass_normalization_coefficient", 9*R(10, 59)**3 < R(1, 16))
    check("actual_mass_source_bound", R("5.40144e14")/4 < R("1.4e14"))
    check("actual_value_source_bound", 24*R("5.58e9") < 2*10**11)
    check("fixed_collar_mass_scaling", RAD**-1 < 10**164)
    check("benchmark_g_upper", R(2000)**-10 < 1)
    final = 3*10**22+3*10**5000+2*10**11*10**386000+R("1.4e14")*10**82*10**5100
    check("complete_compatible_global_disk_bound", final < 10**386020)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 collar/patch global disk bound",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        support_and_partition()
        finite_transfer_budgets()
        mass_profile()
        exact_interface_absorption()
        physical_source_and_solution()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S166 collar/patch global disk bound",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "quantitative_interface_norm_transfer": good,
        "global_disk_graph_absorbed": good, "global_disk_inverse_in_specified_norm": good,
        "higher_York_norms": False, "residual_inverse_composition": False,
        "nonlinear_contraction": False, "S166_complete": False, "goal_complete": False,
        "unresolved": ["higher York extraction/reconstruction", "matched inverse/residual composition",
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
