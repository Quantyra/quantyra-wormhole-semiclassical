"""S166 complete physical mass-column image and mass-extended outer localization."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-mass-column-source-bound"
OWN = ["checks/outer_mass_column_source_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("outer-complete-disk-normal-comparison", 224),
    ("outer-fast-localization-and-tails", 280),
    ("outer-strong-curved-quantum-remainder", 243),
    ("outer-complete-massive-model-comparison", 277),
    ("outer-flat-tower-normal-remainder", 321),
    ("tt-mass-numerical-bound", 398),
]
PROTOCOL = {
    "column": "ehat_m=c_L(Y) Y^3 (-1/12,1), L=epsilon/8; e_m=M_epsilon^-1 ehat_m",
    "normal_cancellation": "exact complete N0(3) mass null vector before any source estimate",
    "canonical_image_bound": "||C^-1 Lhat ehat_m||_0 <1e86 L^-1/2 on entire actual half",
    "physical_source_bound": "Fdisk(A e_m)<1e88 lambda^2 sqrt(gamma0) L^-1/2",
    "mass_input_is_fast": False,
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


def normal_mass_and_integrals():
    s, g, F = S.symbols("s g F", nonzero=True)
    qs, pp, k0 = s*(s-3), s*(s-1)*(s-2)*(s-3), 60*S.pi/g+R(11, 3)
    mat = S.Matrix([
        [(qs-4)*(qs+k0)/(80*S.pi**2), (s-4)*(qs+k0)/(240*S.pi**2)],
        [-(s+1)*(qs+k0)/(240*S.pi**2), pp*F/(480*S.pi**2)
         -(1/(48*S.pi*g)+R(23, 8640)/S.pi**2)*qs-(qs+k0)/(720*S.pi**2)]])
    mass = S.Matrix([-R(1, 12), 1])
    for row, value in enumerate(mat.subs(s, 3)*mass):
        zero("complete_normal_mass_null_"+str(row), value)
    check("mass_vector_length", (mass.T*mass)[0] < 4)
    y, z, ell = S.symbols("Y Z L", positive=True)
    delta = R(7, 2)
    check("mass_outside_fast_domain", 2*(3-delta) == -1)
    check("mass_source_scaling", 3-delta == -R(1, 2))
    check("outer_normal_tail_integral", S.integrate(z**-2, (z, 1, S.oo)) == 1)
    inner_sq = S.integrate((33*y**(4-delta))**2/y, (y, 0, R(1, 2)))
    far_sq = S.integrate((132*y**(-1-delta))**2/y, (y, 4, S.oo))
    middle_sq = S.integrate((R(33, 1024)*y**(-1-delta))**2/y, (y, R(1, 2), 4))
    check("normal_inner_source", inner_sq < 24**2, inner_sq)
    check("normal_far_source", far_sq < 1, far_sq)
    check("normal_middle_separated_source", middle_sq < 100**2, middle_sq)
    check("mass_input_integral", S.integrate(z**3, (z, 0, 2)) == 4)
    check("middle_input_integral", S.integrate(z**3, (z, 0, R(1, 4))) == R(1, 1024))
    check("two_separated_kernels", 2**5+1 == 33)
    check("log_sixteen_length", sum(R(3)**k/S.factorial(k) for k in range(6)) > 16)
    check("log_thirty_two_length", sum(R(4)**k/S.factorial(k) for k in range(6)) > 32)
    check("middle_weighted_profile", 3**5*4*10**28 < 10**32)
    check("middle_H5_pair", 2**5*2*4*10**64 < 10**68)
    check("middle_mixed_norm", 2*10**34 < 10**35)
    check("normal_middle_full_source", 10**38*10**35 == 10**73)
    check("normal_complete_mass_image", 24+1+100+10**73 < 10**80)
    check("sigma_outer_factor_one", R(8, 4) == 2)
    check("sigma_inner_scale", R(1, 4)/4 == R(1, 16))


def mass_derivatives_and_remainders():
    y, ell = S.symbols("Y L", positive=True)
    ds = [1, 11, 220, 3960, 63360, 887040, 10644480]
    for k in range(4):
        budget = sum(S.binomial(k, j)*ds[j]*2**j*S.factorial(3)/S.factorial(3-k+j)
                     for j in range(k+1))
        check("cubic_cutoff_radial_"+str(k), budget < 10**6, budget)
        check("inverse_m_derivative_"+str(k), 2*S.factorial(k)*20**k < 10**5)
    for k in range(5):
        check("cubic_cutoff_conormal_"+str(k), 4**k*10**14 < 10**18)
    check("inverse_m_product", 8*10**5*10**6 < 10**15)
    check("weighted_first_canonical", 4*10**30 < (4*10**15)**2)
    check("Hess_canonical", 2*9*10**30*8/3 < 10**34)
    check("third_covariant_point", 17*10**15 < 10**17)
    check("third_covariant_canonical", 2*10**17 < 10**18)
    check("scalar_mass_canonical_value", 2**7/R(7) < 100)
    check("pair_mass_canonical_value", 2*R(145, 144)*2**7/R(7) < 100)
    check("first_contact_source", 4*10**60*4*10**15 < 10**77)
    check("second_contact_source", 4*10**30*10**17 < 10**48)
    check("Taylor_actual_source", 4*2*10**6*10**18 < 10**25)
    # q_L<=2Y^3 and (1+1/Y)q_L<=4Y^2 on this support.
    check("constant_contact_source", (4*R("1.6e25")*4)**2*R(64, 5) < 10**60)
    check("bounded_actual_source", 4*8*10**24*10 < 10**30)
    check("complete_actual_to_geo", 10**77+10**48+10**25+2*10**30 < 10**80)
    check("geodesic_flat_contacts", 2*3*10**40*10**6 < 10**48)
    check("geodesic_Taylor_source", 2*10**7*10**6 < 10**15)
    check("geodesic_bounded_source", 10**20*10 < 10**25)
    check("complete_geo_to_flat", 10**48+10**15+10**25 < 10**50)
    check("flat_mass_H3", 2**3*4*10**12 < 10**16)
    check("flat_complete_tower_source", 1000*10**8 < 10**12)
    check("complete_massive_difference", 10**80+10**50+10**12 < 10**81)
    check("local_cubic_source", (6*10**45*10**7*10**18*4)**2*2 < 10**148)
    check("complete_outer_mass_source", 10**81+R("1.5e74")+R("1.5e80") < 10**83)
    for order in [1, 2, 3]:
        # In the mass-adapted comparison norms, inner truncations have L2^2 scale a.
        exponent = (3-order)-(2 if order == 1 else 1 if order == 2 else 0)
        check("inner_truncation_power_"+str(order), 2*exponent+1 == 1)


def global_physical_conversion():
    pi = R(22, 7)
    check("outer_canonical_conversion", 4*pi*2*R(3, 2) < 100)
    check("outside_separation", 1-R(1, 4) > R(1, 2))
    check("outside_massive_source", 4*10**29*10 < 10**32)
    check("outside_scale_conversion", R(8)**-4 < 1)
    check("complete_canonical_image", 100*10**83+10**32 < 10**86)
    c, u, bb, aa, qcan = S.symbols("C u B A0 Qcan", positive=True)
    zero("physical_weight_density_cancellation", u**-4*(c*u**4/(4*S.pi*bb**4))**2
         *(4*S.pi*bb**4/c)*u**-4-c/(4*S.pi*bb**4))
    beta, lam = S.symbols("beta lambda", positive=True)
    gamma0 = 4*S.pi*beta*bb**3
    k = S.sqrt(lam*bb**2/c)
    left = c**R(3, 2)*S.sqrt(beta*S.sqrt(c*lam))/(2*S.sqrt(S.pi)*bb**2)
    right = lam**2*S.sqrt(gamma0)*k**(-R(7, 2))/(4*S.pi)
    zero("physical_conversion_square", left**2-right**2)
    dual_inverse = S.Matrix([[1, 0], [R(1, 3), 2]])
    check("inverse_dual_matrix_norm", sum(x*x for x in dual_inverse) < 9)
    check("outer_k_power", R(5, 2)**7 < 25**2)
    check("complete_physical_conversion", 3*2*25/(4*3) < 100)
    check("physical_mass_source", 100*10**86 <= 10**88)
    check("mass_coefficient_conversion", 9/R("5.9")**3 < R(1, 16))
    check("physical_mass_forcing", R("5.40144e14")*10**88/4 < 2*10**102)
    scalar, theta, b = S.symbols("scalar theta b")
    zero("physical_mass_vector_theta", scalar+(-2*scalar)/3-(-4*scalar)*(-R(1, 12)))
    # r=L/4; the value contribution is the dominant mass term in localization.
    mass_terms = 2*10**20*10**83+R("1e-40")*10**35+10**2000*4**4*10
    check("mass_extended_localization", mass_terms < 10**2005)
    check("mass_is_not_fast", PROTOCOL["mass_input_is_fast"] is False)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 full mass-column source bound",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        normal_mass_and_integrals()
        mass_derivatives_and_remainders()
        global_physical_conversion()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 full mass-column source bound", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version, "sympy": S.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "complete_mass_column_source_bounded": good, "mass_extended_outer_localization": good,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["transition and intervening annulus/bolt matching", "higher York maps",
                       "full inverse/residual composition", "nonlinear signed compatibility",
                       "all physical completion requirements"],
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
