"""S166 higher genuine York maps in the fixed scalar fourth/log atlas norms.

Symbolic identities and finite reserves support the accompanying ordinary
global proof; they do not formally certify its functional-analytic steps.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
from sympy.functions.combinatorial.numbers import stirling

ROOT = Path(__file__).resolve().parents[1]
STEM = "quantitative-higher-york-isomorphism"
OWN = ["checks/quantitative_higher_york_isomorphism.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("collar-patch-global-disk-bound", 281),
    ("weighted-york-source-conversion", 293),
    ("tt-mass-numerical-bound", 398),
    ("outer-mass-column-source-bound", 244),
    ("outer-strong-curved-quantum-remainder", 243),
    ("stationary-weighted-tt-inverse", 52),
    ("stationary-global-tt-reduction", 44),
]
R = S.Rational
EPS, TINY = R(1, 10**160), R(1, 10**170)
RAD, ELL = EPS/32, EPS/8
ROWS = []
PROTOCOL = {
    "scope": "genuine linear seed York maps in explicit scalar fourth/log norms",
    "fast_maps": "both directions below 1e24000",
    "scalar_mass_column": "c_L(Y)u(Y)^3; m_Y=-3k^3 m_c",
    "mass_fast_discrepancy": "below 1e13000 L^-1/2",
    "forward": "Ynorm<=1e25000 G",
    "reverse": "G<=1e40000 Ynorm",
    "compatible_inverse": "Ynorm(Rcal^-1 A_m^-1 f)<=1e412000 lambda^-2 gamma0^-1/2 Fdisk",
    "residual_inverse_composition": False,
    "nonlinear_candidate_regularity": False,
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


def weighted_energy():
    aa, rr, vv, zp, zz, tt = S.symbols("A R V Zprime Z T")
    wp, rp = aa+2*rr, aa*rr-rr**2
    # Integration of 2R Z Zprime gives -(Rprime+wR)Z^2.
    conjugate_potential = vv-6*rr**2+3*rp+3*wp*rr
    zero("weighted_conjugated_energy", conjugate_potential-3*(rp+wp*rr)-(vv-6*rr**2))
    square_extra = R(3, 2)*(R(9, 4)*rr**2-R(3, 2)*(rp+wp*rr))
    zero("weighted_square_completion", square_extra+R(9, 2)*aa*rr-R(57, 8)*rr**2+6*rr**2)
    # -W(Tprime+4RT), integrated against Z, with Wprime/W=2R.
    zero("weighted_source_integration", tt*(zp+(2*rr+wp-4*rr)*zz)-tt*(zp+aa*zz))
    check("projection_ratio_from_U", (R(399, 100)-2)/2 == R(199, 200))
    check("weighted_positive_R_square", R(9, 2)*R(199, 200)-R(33, 8) == R(141, 400))
    check("rounded_R_square", R(141, 400) > R(7, 20))
    cost = R(2, 3)+R(2, 3)+R(9, 4)/R(7, 20)
    check("energy_source_Cauchy", cost == R(163, 21) and cost < 8)
    check("projection_vector_low", 8/R("1.49") < 3**2)
    check("projection_A_low", 8/R("1.5") < 3**2)
    check("projection_R_low", 8/R(".35") < 5**2)
    check("projection_derivative_low", 8*(R(2, 3)+R(49, 4)/R(".35")) < 17**2)
    check("extraction_low_system", 3+17+3+2*5+R(2, 3) < 500)
    check("return_vector_low", R(2)*11**2 < 16**2)
    check("return_derivative_low", 11+11 == 22)
    check("return_R_low", R("3.2")*11 == R("35.2"))
    check("return_low_system", 16+22+11+4*R("35.2")+1 < 500)
    # Arbitrary weighted vector forcing: E<=||source|| ||Z||.
    check("forced_vector_low", 1/R("1.49") < 1)
    check("forced_square_and_A_low", 1/(R("1.49")*R("1.5")) < 1)
    check("forced_R_low", 1/(R("1.49")*R(".35")) < 2**2)
    check("forced_divergence_low", 1+1+1+R(3, 2)*2 < 500)


def geometry_jets():
    u, nu, alpha, beta, ee = S.symbols("u nu alpha beta e")
    ff = 1+nu*u**2-2*alpha*u**3+beta*u**4-2*ee*u**5+ee*u**6
    aa = ff-u*S.diff(ff, u)/2
    kk1 = -ff+u*S.diff(ff, u)-u**2*S.diff(ff, u, 2)/2
    # Derive the proper warped-curvature coefficients from F=lambda r^2 f(B/r).
    radius, lam, bb = S.symbols("r lambda B", positive=True)
    fr = lam*radius**2*ff.subs(u, bb/radius)
    zero("geometric_a_polynomial", (S.diff(fr, radius)/(2*lam*radius)).subs(radius, bb/u)-aa)
    zero("geometric_K1_polynomial", (-S.diff(fr, radius, 2)/(2*lam)).subs(radius, bb/u)-kk1)
    bounds = {nu: R(".001"), alpha: R("2.002"), beta: R("3.002"), ee: R("1e-15")}
    for name, poly in [("f", ff), ("a", aa), ("K1", kk1), ("K2", -aa), ("K3", nu*u**2-ff)]:
        for n in range(7):
            total = sum(sum(abs(c)*S.prod(bounds[v]**power for v, power in zip(bounds, mon))
                            for mon, c in S.Poly(coef, *bounds).terms())
                        for coef in S.Poly(poly, u).all_coeffs())
            check(f"polynomial_{name}_jet_{n}", total < 10**9)
            poly = S.expand(u*S.diff(poly, u))
    b0, k0 = 10**9, 10**10
    for n in range(7):
        check("tensor_induction_step_"+str(n), 3*(n+1)*b0**2*k0**n*S.factorial(n)
              <= b0*k0**(n+1)*S.factorial(n+1))
        raw = b0*k0**n*S.factorial(n)
        curv = 16*S.binomial(n+2, 2)*b0**3*k0**n*S.factorial(n)
        check("geometric_global_jet_"+str(n), max(raw, curv) < 10**110)
        check("collar_geometric_jet_reserve_"+str(n), raw < S.factorial(n)*10**(60*(n+1)))
    check("quadratic_system_coefficient", 12*2**5*10**220 < 10**230)


def first_order_systems():
    dd, yp, yy, aa, rr, tt, cc, ss, ww, ric, tp, cp = S.symbols(
        "D Yprime Y A R T c S W Ric Tprime cprime")
    div = yp+(aa+2*rr)*yy
    theta = ss+tt/3
    tau = ss+tt/2-dd/4
    cextract = -tt+dd/2-2*rr*yy
    psi = ww*(dd-2*tt/3)
    zz, uu = ww*yy, ww*tt
    zero("extract_scalar_output", ww*(tau-theta)+psi/4)
    zero("extract_shape_output", ww*cextract-(-2*uu/3+psi/2-2*rr*zz))
    zero("extract_divergence", ww*dd+2*rr*zz-(2*uu/3+psi+2*rr*zz))
    dp = R(2, 3)*(tp+4*rr*tt-2*ric*yy)
    psip = 2*rr*psi+ww*(dp-2*tp/3)
    zero("extract_gradient_equation", psip-(2*rr*psi+R(8, 3)*rr*uu-R(4, 3)*ric*zz))
    # The return has eta=-Y, b=-Yprime+AY-c, so use its original formulas.
    br = -yp+aa*yy-cc
    sc = cc/2+rr*yy
    tc = -(br+3*cc)/2+(aa-rr)*yy
    zero("return_T_kinematics", tc-(div/2-2*rr*yy-cc))
    psir = ww*(div+2*rr*yy+cc)
    zero("return_theta_output", ww*(sc+tc/3)-psir/6)
    zero("return_T_output", ww*tc-(-R(3, 2)*ww*cc+psir/2-3*rr*zz))
    zero("return_divergence", ww*div+2*rr*zz-(psir-ww*cc))
    rp = aa*rr-rr**2
    # Solve the original projection after T=D/2-2RY-c.
    dpreturn = -2*rp*yy-2*rr*yp-cp+2*rr*dd-8*rr**2*yy-4*rr*cc-2*ric*yy
    psi_symbolic = ww*(dd+2*rr*yy+cc)
    psiprime = 2*rr*psi_symbolic+ww*(dpreturn+2*rp*yy+2*rr*yp+cp)
    zero("return_gradient_equation", psiprime-(4*rr*psi_symbolic-6*rr*ww*cc
         -2*ric*zz-12*rr**2*zz))
    forcing = S.Symbol("forcing")
    dpforce = R(2, 3)*(forcing-2*ric*yy)
    zero("forced_gradient_equation", 2*rr*ww*dd+ww*dpforce
         -(2*rr*ww*dd+R(2, 3)*ww*forcing-R(4, 3)*ric*zz))


def hodge_and_induction():
    kk = 10**110
    for j in range(1, 6):
        comm = 2**(j-1)+16*(j+1)*(j-1)*2**max(j-2, 0)+16*(j+1)
        check("Hodge_commutator_reserve_"+str(j), comm*(1+kk) < 10**120)
        check("second_derivative_energy_reserve_"+str(j), 24*(j+1)*kk
              < (8*(j+1)*(1+kk))**2)
    for j in range(6):
        # Both gradient and divergence identities, product terms and Hodge remainder.
        budget = 20*(j+2)*2**j*10**230+10**120+1
        check("first_order_induction_budget_"+str(j), budget < 10**235)
        check("higher_induction_step_"+str(j), 10**235*(10**(240*j+3)+1)
              < 10**(240*(j+1)+3))
    check("low_induction_start", 500 < 10**3)
    for order in [3, 5]:
        budget = (order+1)*2**order*10**230*(10**(240*order+3)+1)
        check("covariant_output_H"+str(order), 3*budget < 10**2000)
    check("forced_mass_correction_H5", 6*2**5*10**110*10**1203 < 10**1400)


def atlas_transfer():
    check("interior_u_lower", RAD/24 > R(1, 10**164))
    check("interior_conformal_upper", R(1, 2)*10**164 < 10**164)
    check("derivative_conversion_base", 10**164*10**170*10**550*10**110 < 10**1000)
    check("weighted_patch_measure_root", R(1, 2)*10**(8*164) < 10**2000)
    budgets = []
    for j in range(6):
        # Bell sums and product of two derivative coefficient families.
        bell = 1 if j == 0 else sum(S.bell(j, m, tuple(S.factorial(n)
                    for n in range(1, j-m+2))) for m in range(1, j+1))
        product = sum(S.binomial(j, n)*S.factorial(n)*S.factorial(j-n) for n in range(j+1))
        check("atlas_Bell_product_"+str(j), bell*product <= S.factorial(j)*100**j)
        budgets.append(10**1000*S.factorial(j)*(100*10**1000)**j)
    for order in [3, 5]:
        endpoint = 10**1000*4*1000*3**order*max(budgets[:order+1])
        check("atlas_both_endpoints_H"+str(order), endpoint < 10**10000)
    k, a0, bb, cc, u, yy, ws, jj, beta = S.symbols("k A0 B C u Y w j beta", positive=True)
    gamma = 4*S.pi*a0*bb**4/(k*cc)
    zero("weighted_physical_patch_measure", (4*S.pi*bb**4/cc)*u**-8/gamma-(k/a0)*u**-8)
    zero("weighted_outer_density", k*(k*yy*ws)**-8*jj*yy*yy**7-k**-7*ws**-8*jj)
    part, ff = S.symbols("partition F", nonzero=True)
    zero("analysis_extension_identity", (1-part)*ff+part*ff-ff)
    period = S.Symbol("period", positive=True)
    zero("actual_time_average", S.integrate(ff, (yy, 0, period))/period-ff)
    n = S.Symbol("n", integer=True, positive=True)
    schur = 32*(1+2*S.summation((n+1)/2**n, (n, 1, S.oo)))
    check("dyadic_log_transfer", 128*schur == 28672)
    check("dyadic_plain_transfer", 512*(1+2*S.summation(2**-n, (n, 1, S.oo))) == 1536)
    check("ambient_endpoint_composition", 10000+2000+10000 == 22000)
    check("logarithmic_actual_maps", 28672*10**22000 < 10**23000)
    check("mixed_fast_map_reserve", 20*10**23000 < 10**24000)


def mass_comparison():
    yy = S.Symbol("Y", positive=True)
    te, ws, aa, jj = (S.Function(name)(yy) for name in ["T", "w", "A", "j"])
    dd = lambda expr: yy*S.diff(expr, yy)
    radial = lambda expr: -ws*dd(expr)
    vv = S.Symbol("V")
    qe = -te/6
    p_qe = R(3, 2)*(radial(radial(qe))+(aa+2*jj)*radial(qe))-vv*qe
    actual = radial(te)+4*jj*te-p_qe
    model = ws**2*dd(dd(te))/4+(ws*dd(ws)/4-(1+(aa+2*jj)/4)*ws)*dd(te)+(4*jj-vv/6)*te
    zero("mass_comparison_operator", actual-model)
    s = S.Symbol("s")
    zero("mass_normal_factorization", s**2/4-R(7, 4)*s+3-(s-3)*(s-4)/4)
    zero("mass_normal_cubic_cancellation", (s**2/4-R(7, 4)*s+3).subs(s, 3))
    te0 = yy**3/2
    div0 = -dd(-te0/6)+3*(-te0/6)
    zero("mass_model_no_scalar_cubic", div0)
    zero("mass_model_shape_cubic", -te0+div0/2-2*(-te0/6)+yy**3/3)
    k, w, im = S.symbols("k w inverse_m", positive=True)
    zero("mass_local_fast_subtraction", -2*(im*yy**3/2)/3+(k*yy*w)**3/(3*k**3)
         -yy**3*(w**3-im)/3)
    check("mass_source_smaller_analytic_domain", 2*ELL < R(".001"))
    check("Cauchy_first_jet", R(".001")/R(".05") == R(".02"))
    check("Cauchy_second_jet", 2*R(".001")/R(".05")**2 == R(".8"))
    ab = R("1.001")+R("1.001")*R(".02")**2/R(".999")
    ap = R(".02")+R(".02")**3/R(".999") \
        +R("1.001")*(R(".02")+R(".02")*R(".8"))/R(".999") \
        +R("1.001")*R(".02")*R(".02")**2/R(".999")**2
    check("analytic_A_modulus", ab < R("1.002"))
    check("analytic_A_derivative", ap < R(".1"))
    vbound = R("1.001")*R(".02")*(R(".1")+2*R(".02"))/2+2*ab**2+4*R("1.001")**2
    check("analytic_V_modulus", vbound < 10)
    # Cauchy coefficient bound on radius .01 and Y<=.001, after two vanishing terms.
    n = S.Symbol("n", integer=True, nonnegative=True)
    series = 10*100**2*S.summation((n+2)**6/R(10)**n, (n, 0, S.oo))
    check("mass_two_power_analytic_jets", series < 10**30)
    z = S.Symbol("z")
    cutoff = 1-(462*z**6-1980*z**7+3465*z**8-3080*z**9+1386*z**10-252*z**11)
    derivative_bounds = [1, 11, 220, 3960, 63360, 887040, 10644480]
    for order in range(7):
        poly = S.Poly(S.diff(cutoff, z, order), z)
        degree = poly.degree()
        bernstein = [sum(poly.nth(j)*S.binomial(i, j)/S.binomial(degree, j)
                         for j in range(i+1)) for i in range(degree+1)]
        check("beta_cutoff_derivative_"+str(order), max(abs(v) for v in bernstein) <= derivative_bounds[order])
        logbound = sum(stirling(order, j, kind=2)*2**j*derivative_bounds[j] for j in range(order+1))
        check("beta_cutoff_log_jet_"+str(order), logbound < 10**14)
        if order < 6:
            check("beta_cutoff_endpoint_regular_"+str(order),
                  S.diff(cutoff, z, order).subs(z, 1) == 0
                  and S.diff(cutoff-1, z, order).subs(z, 0) == 0)
    for order in range(7):
        cubic = sum(S.binomial(order, j)*10**14*3**(order-j) for j in range(order+1))
        product = 2**order*10**30*cubic
        check("mass_cubic_with_inverse_density_"+str(order), product < 10**60)
    for order in range(6):
        # Output error coefficients (including two-power cancellation) times up to six input jets.
        budget = 10*2**order*10**30*10**60
        check("mass_scalar_radial_error_"+str(order), budget < 10**100)
        if order <= 4:
            check("mass_vector_source_radial_"+str(order), budget < 10**100)
        bell = 1 if order == 0 else sum(S.bell(order, j, tuple(S.factorial(l)
                for l in range(1, order-j+2))) for j in range(1, order+1))
        check("weighted_covariant_chain_"+str(order), 2**order*S.factorial(order)**2*bell*10**600 < 10**700)
    check("collar_weighted_measure", R(".4")**-7*R(".999")**-8*R("1.001") < 10**5)
    ell = S.Symbol("L", positive=True)
    zero("mass_inner_square_integral", S.integrate(yy**2, (yy, 0, ell))-ell**3/3)
    zero("mass_transition_square_integral", S.integrate(yy**-2, (yy, ell, 2*ell))-1/(2*ell))
    check("mass_covariant_square_budget", 6*2**5*10**5*(10**100*10**700)**2*(ELL**4/3+R(1, 2)) < 10**2000)
    # Include both local errors, forced correction, atlas conversion, log<=2H5.
    check("mass_full_fast_discrepancy", 4*10**10000*(10**1400*10**1000+10**1000) < 10**13000)
    momentum = S.Symbol("p", nonnegative=True)
    zero("log_by_linear_derivative", S.diff(1+momentum-S.log(2+momentum), momentum)
         -(1+momentum)/(2+momentum))
    check("log_by_linear_at_zero", 1+1+R(1, 2) > 2,
          "The first three terms of exp(1) exceed 2, hence log(2)<1.")
    zero("H4log_by_H5", 2*(1+momentum**2)-(1+momentum)**2-(momentum-1)**2)


def final_assembly():
    check("cutoff_L_inverse_root", ELL**-1 < 10**164)
    check("mass_extraction_factor", 1/(3*R(".4")**3) < 6)
    check("mass_return_factor", 3*R(".5")**3 < 1)
    forward = 10**24000+10**13000*10**82+6
    reverse = 10**24000+10**24000*10**13000*10**82+1
    check("higher_York_forward_complete", forward < 10**25000)
    check("higher_York_reverse_complete", reverse < 10**40000)
    check("compatible_TT_scalar_inverse", 386020+25000 < 412000)
    # Linear mass identity independent of a choice of fast norm.
    k, e, tau, cf, cm = S.symbols("k e tau fastc massc")
    zero("mass_reverse_discrepancy_identity", (cf-cm/(3*k**3))+cm/(3*k**3)-cf)
    check("normal_mass_not_fast", 2*(3-R(7, 2)) == -1)
    check("projection_return_tail_fast", 2*(4-R(7, 2)) > 0)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 quantitative higher York isomorphism",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        weighted_energy()
        geometry_jets()
        first_order_systems()
        hodge_and_induction()
        atlas_transfer()
        mass_comparison()
        final_assembly()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S166 quantitative higher York isomorphism",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "higher_York_norms": good, "explicit_mass_representative_comparison": good,
        "global_compatible_TT_scalar_bound": good,
        "residual_inverse_composition": False, "nonlinear_contraction": False,
        "S166_complete": False, "goal_complete": False,
        "unresolved": ["matched inverse/residual composition", "S166 completion audit",
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
