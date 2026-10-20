"""S166 actual residual in the physical forcing norm and compatible composition."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "matched-residual-inverse-composition"
OWN = ["checks/matched_residual_inverse_composition.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("quantitative-higher-york-isomorphism", 390),
    ("comparison-path-complete-equation-residual", 329),
    ("comparison-path-complete-source-response", 305),
    ("general-background-newton-operator", 174),
    ("stationary-scalar-outgoing-inverse", 107),
]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual comparison forcing and full seed-compatible inverse composition",
    "forcing": "lambda^-2 gamma0^-1/2 Fdisk(f_ell)<2e-8 ell^2/g<2.2e25 ell^2",
    "projection": "orthogonal projection in the exact forcing Hilbert space, norm one",
    "disk_correction": "G(h)<2*10^386012 ell^2/g<3*10^386045 ell^2",
    "York_correction": "Ynorm(Rcal^-1 h)<2*10^411992 ell^2/g<3*10^412025 ell^2",
    "candidate_Hessian_inverse": False,
    "signed_compatibility_evaluated": False,
    "useful_Newton_certificate": False,
    "fixture_digits": 80, "fixture_h": ["1e-5", "1e-30"],
    "fixture_t": ["1+2h", "1.1", "4", "50"], "fixture_ell": ["0", ".4", "1"],
    "fixture_amplitude": ".02", "fixture_normalization": ".02^2 min(t-1,1)",
    "fixture_tolerance": "1e-30(1+abs(target))", "fixtures_are_physical_states": False,
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


def exact_orbit_algebra():
    xx, xxpp, sig, potential = S.symbols("X Xpp sigma potential", nonzero=True)
    euler_sigma = -2*xxpp+2*S.exp(2*sig)*potential
    zero("reduced_orbit_variation", -euler_sigma/(2*S.exp(2*sig)*xx)
         -(S.exp(-2*sig)*xxpp-potential)/xx)
    t = S.Symbol("t", positive=True)
    ff = S.Function("f")(t)
    hh = t**2*S.diff(t*ff, t)
    delta = lambda fn: t**2*(ff*S.diff(fn, t, 2)+S.diff(ff, t)*S.diff(fn, t))
    zero("seed_orbit_identity", delta(t**2)-2*hh)
    tau = 1-1/t
    drift = delta(tau)+4*t*ff*S.diff(tau, t)
    zero("orbit_drift_identity", drift-S.diff(ff, t)-2*ff/t)
    ell, u, nu, tt, kk, h0, h1, hn, star, lnu, ltau, gv2, gt2, vv, aa = S.symbols(
        "ell u nu T k H Hprime Hnew Ostar Lnu Ltau gradV2 gradT2 V a")
    bb = kk-1-ell*vv*kk
    lt = ell*kk*lnu+bb*ltau+tau*ell**2*kk*gv2
    original = aa*(2*h0*tt**2+2*t**2*tt*lt+2*t**2*gt2)-2*hn
    tangent = 4*(nu-u)*h0+2*t**2*lnu-2*t*h1*nu
    exact = ell*(aa*tt*kk-1)*star \
        +2*h0*(aa*(tt**2+2*ell*tt*kk*(u-nu))-1) \
        +2*ell*aa*tt*kk*t*h1*nu-2*(hn-h0) \
        +2*aa*t**2*tt*bb*ltau+2*aa*t**2*tt*tau*ell**2*kk*gv2+2*aa*t**2*gt2
    zero("complete_finite_orbit_decomposition", (original-ell*tangent-exact).subs(star, tangent))
    # Full first-order drift product rule with nu=tau V and T=1+tau(k-1).
    tau0, tau1, tau2, v1, v2, drift0 = S.symbols("tau tau1 tau2 V1 V2 drift")
    kt = ell*kk*v1
    ktt = kk*(ell*v2+ell**2*v1**2)
    direct = tau2*(kk-1)+2*tau1*kt+tau0*ktt+drift0*(tau1*(kk-1)+tau0*kt)
    lnu_value = tau2*vv+2*tau1*v1+tau0*v2+drift0*(tau1*vv+tau0*v1)
    zero("radius_exponential_drift_product", direct-(ell*kk*lnu_value
         +bb*(tau2+drift0*tau1)+tau0*ell**2*kk*v1**2))
    xrat, erat, drat, arat = S.symbols("xrat erat drat arat")
    fpoly = erat+xrat*t**2-2*arat/t+(erat+3*xrat+drat)/t**2-2*drat/t**3+drat/t**4
    hp = 3*xrat*t**4+erat*t**2-(erat+3*xrat+drat)+4*drat/t-3*drat/t**2
    zero("actual_H_polynomial", t**2*S.diff(t*fpoly, t)-hp)
    check("H_segment_bound", R("1.001")*(3/R("5.9")+R(".001")+R(".501")+7*R("1e-24")) < R("1.1"))
    check("H_first_segment_bound", R("1.001")*(12/R("5.9")+R(".002")+10*R("1e-24")) < R("2.1"))
    check("H_second_segment_bound", R("1.001")*(36/R("5.9")+R(".002")+26*R("1e-24")) < R("6.3"))
    check("drift_cold", 5+3*R("1.01")*R(".001") < 6)
    check("drift_middle", 3+3 == 6)
    check("drift_outer", 8+15/R("5.9") < 11)


def source_and_regional_bounds():
    source_rows = [
        ("cold", R("1.7e-6"), R("1e-14"), R("2.6e-10"), R("3.2e-16"), R("3.401e-6"), R("5.201e-10")),
        ("middle", R("2.6e-9"), R("1e-15"), R("2.2e-12"), R("3.2e-17"), R("5.201e-9"), R("4.401e-12")),
        ("outer", R("4.7e-7"), R("1e-12"), R("4.9e-11"), R("3.2e-14"), R("9.401e-7"), R("9.807e-11")),
    ]
    for name, av, tv, ad, td, ob, od in source_rows:
        check("source_orbit_seed_"+name, 2*(av+tv) < ob)
        check("source_orbit_first_"+name, 2*(ad+td) < od)
    check("outer_Ostar_enlarges_middle", R("5.201e-9") < R("9.401e-7"))
    tm, am = R("1.001"), R(".00331")
    uc, nc, w, b = R("4e-6"), R("1.2e-8"), R("1.2e-8"), R("3e-7")
    fac = 2*uc*am+w+R("1.001")*nc*am
    cold = R("1.01")*fac*R("3.401e-6") \
        +2*R("1.1")*tm**4*(10*(uc+nc)**2*am+5*(uc+nc)*w+2*w**2) \
        +R("2.02")*R("2.1")*tm**4*nc*fac \
        +R("1.002")*R("2.1")*tm**4*w**2 \
        +R("1.004")*R("6.3")*tm**4*nc**2*am \
        +R("1.01")*w**2*6 \
        +R("2.02")*R("1.01")*tm*b**2 \
        +R("2.02")*R("1.01")*R(".001")*(w+b)**2
    check("geometric_orbit_cold", cold < R("3e-12"), S.N(cold, 14))
    t, ss = R(4), R(3)
    u = n = w = R("1.2e-8")
    b = R("3e-7")
    fac = 2*u+w+R("1.001")*n
    middle = R("1.01")*fac*R("5.201e-9") \
        +2*R("1.1")*t**4*(10*(u+n)**2+5*(u+n)*w+2*w**2) \
        +R("2.02")*R("2.1")*t**4*n*fac \
        +R("1.002")*R("2.1")*t**4*w**2 \
        +R("1.004")*R("6.3")*t**4*n**2 \
        +R("1.01")*w**2*6*t**3 \
        +R("2.02")*t**3*ss*b**2 \
        +R("2.02")*ss**2*(w+t*b)**2
    check("geometric_orbit_middle", middle < R("7e-11"), S.N(middle, 14))
    u, n, w, b = R("8e-9"), R("2e-9"), R("2e-9"), R("1.4e-8")
    fac = 2*u+w+R("1.001")*n
    outer = R("1.01")*fac*R("9.401e-7")/16 \
        +2*R("1.1")*(10*(u+n)**2+5*(u+n)*w+2*w**2) \
        +R("2.02")*R("2.1")*n*fac+R("1.002")*R("2.1")*w**2 \
        +R("1.004")*R("6.3")*n**2+R("1.01")*w**2*11/4 \
        +R("2.02")*5/R("5.9")*R(4, 3)*b**2 \
        +R("2.02")*5/R("5.9")*(w/4+b)**2
    check("geometric_orbit_outer", outer < R("5e-15"), S.N(outer, 14))
    no, ns, ell, dt, da = S.symbols("Norb Nangular ell DeltaTrace DeltaAngular")
    rv, rs = no-2*ell*(dt-da), ns-ell*da
    zero("complete_trace_cancellation", rv+2*rs-(no+2*ns-2*ell*dt))
    for name, no, ns, dt, dv, phib, vb in [
        ("cold", "3e-12", "2e-12", "3.2e-16", "5.201e-10", "7.1e-12", "5.3e-10"),
        ("middle", "7e-11", "6e-11", "3.2e-17", "4.401e-12", "2e-10", "8e-11"),
        ("outer", "5e-15", "2e-15", "3.2e-14", "9.807e-11", "7.4e-14", "1e-10"),
    ]:
        check("complete_trace_"+name, R(no)+2*R(ns)+2*R(dt) < R(phib))
        check("complete_orbit_"+name, R(no)+R(dv) < R(vb))
    check("outer_trace_source_covers_middle", R("3.2e-17") < R("3.2e-14"))
    check("outer_orbit_source_covers_middle", R("4.401e-12") < R("9.807e-11"))
    check("global_trace_cold_envelope", R("7.1e-12")*am < R("2e-10"))
    check("global_orbit_cold_envelope", R("5.3e-10")*am < R("1e-10"))
    check("global_trace_outer_envelope", R("7.4e-14") < R("2e-10"))
    check("global_orbit_middle_envelope", R("8e-11") < R("1e-10"))


def forcing_norm_and_composition():
    t, tt, ell, u, c, bb, gg, lam, gamma = S.symbols("t T ell u C B G lambda gamma0", positive=True)
    ratio = tt**4*S.exp(2*ell*u)
    zero("exact_density_radius_cancellation", ratio*(t*tt)**-4-t**-4*S.exp(2*ell*u))
    zero("weighted_source_measure", (t**2*t**-4)**2*t**2-t**-2)
    check("volume_exponential_bound", 1/(1-2*R("1.2e-8")) < R("1.001"))
    moment = sum(S.binomial(8, j)*S.factorial(j) for j in range(9))
    check("eighth_log_moment", moment == 109601)
    xx = S.Symbol("x", nonnegative=True)
    zero("log_moment_integral", S.integrate(S.exp(-xx)*(1+xx)**8, (xx, 0, S.oo))-moment)
    check("forcing_integrated_coefficient", R("1.001")**2*5*R("1e-20")*moment < R("7.5e-8")**2)
    zero("physical_inverse_units", c/(gg*bb**2*lam**2)-(c/(lam*bb**2))/(gg*lam))
    check("normalized_physical_residual", R("7.5e-8")*R("6.1")/24 < R("2e-8"))
    check("benchmark_inverse_g", R(2000)**10/R(".999") < R("1.1e33"))
    check("benchmark_residual_constant", R("2e-8")*R("1.1e33") == R("2.2e25"))
    check("disk_composition_exponent", 386020-8 == 386012)
    check("York_composition_exponent", 412000-8 == 411992)
    check("disk_fixed_benchmark_composition", R("2.2")*10**(386012+33) < 3*10**386045)
    check("York_fixed_benchmark_composition", R("2.2")*10**(411992+33) < 3*10**412025)
    check("current_majorant_not_a_small_step_certificate", 3*10**412025 > 1,
          "This is a limitation of the stated upper bound, not a lower bound on the actual step.")


def weighted_projection_fixture():
    metric = S.diag(2, 3, 5)
    chi = S.Matrix([1, -2, 3])
    riesz = metric.inv()*chi
    norm2 = (chi.T*riesz)[0]
    proj = S.eye(3)-riesz*chi.T/norm2
    basis = S.Matrix([[2, -3], [1, 0], [0, 1]])
    check("weighted_adjoint_nonzero", norm2 > 0 and basis.rank() == 2)
    for name, mat in [
        ("projection_idempotent", proj*proj-proj),
        ("projection_selfadjoint", proj.T*metric-metric*proj),
        ("projection_actual_range", chi.T*proj),
        ("projection_norm_one_decomposition", metric-proj.T*metric*proj-chi*chi.T/norm2),
        ("weighted_range_inverse", basis*(basis.T*metric*basis).inv()*basis.T*metric-proj),
    ]:
        check(name, mat == S.zeros(*mat.shape))
    forcing = S.Matrix([7, -5, 11])
    correction = -(basis.T*metric*basis).inv()*basis.T*metric*forcing
    residual = forcing+basis*correction
    check("projected_linear_residual_exact", residual == riesz*(chi.T*forcing)[0]/norm2)
    zero("projection_Pythagoras", (forcing.T*metric*forcing)[0]
         -((proj*forcing).T*metric*(proj*forcing))[0]-(chi.T*forcing)[0]**2/norm2)
    check("nonzero_obstruction_retained", (chi.T*forcing)[0] != 0)


def manufactured_orbit_fixtures():
    mp.mp.dps = 80
    ee, dd = mp.mpf('.001'), mp.mpf('.00002')
    xr = (1-ee-dd)/6
    amp = mp.mpf('.02')
    fu = lambda t: amp*(t-1)/(1+t)**3
    fv = lambda t: amp/(1+t)**3
    hp = lambda t: 3*xr*t**4+ee*t**2-(ee+3*xr+dd)+4*dd/t-3*dd/t**2
    fext = lambda t: (t-1)**2/t**2*(ee+xr*(t**2+2*t+3)+dd/t**2)
    for hi, hh in enumerate([mp.mpf('1e-5'), mp.mpf('1e-30')]):
        root = 1+hh
        thermal = root*fext(root)
        f0 = lambda t: fext(t)-thermal/t
        delta = lambda fn, t: t**2*(f0(t)*mp.diff(fn, t, 2)+mp.diff(f0, t)*mp.diff(fn, t))
        for ti, t in enumerate([1+2*hh, mp.mpf('1.1'), mp.mpf(4), mp.mpf(50)]):
            tau = (t-1)/t
            nu_fn = lambda z: (z-1)/z*fv(z)
            nu = nu_fn(t)
            ln = delta(nu_fn, t)+4*t*f0(t)*mp.diff(nu_fn, t)
            star = 4*(nu-fu(t))*hp(t)+2*t**2*ln-2*t*mp.diff(hp, t)*nu
            for ei, ell in enumerate([mp.mpf(0), mp.mpf('.4'), mp.mpf(1)]):
                ratio = lambda z: 1+(z-1)/z*mp.expm1(ell*fv(z))
                radius = lambda z: z*ratio(z)
                product = mp.exp(-2*ell*fu(t))*delta(lambda z: radius(z)**2, t)-2*hp(radius(t))
                freal = lambda z: f0(z)*mp.diff(radius, z)**2/(ratio(z)**2*mp.exp(2*ell*fu(z)))
                lapse = lambda z: 2*mp.log(ratio(z))+2*ell*fu(z)-mp.log(mp.diff(radius, z))
                areal = 2*radius(t)**2*(freal(t)+radius(t)*(mp.diff(freal, t)
                         +freal(t)*mp.diff(lapse, t))/mp.diff(radius, t))-2*hp(radius(t))
                a, k, tt = mp.exp(-2*ell*fu(t)), mp.exp(ell*fv(t)), ratio(t)
                br = k-1-ell*fv(t)*k
                lt = mp.diff(f0, t)+2*f0(t)/t
                gv2, gt2 = t**2*f0(t)*mp.diff(fv, t)**2, t**2*f0(t)*mp.diff(ratio, t)**2
                decomposition = ell*(a*tt*k-1)*star \
                    +2*hp(t)*(a*(tt**2+2*ell*tt*k*(fu(t)-nu))-1) \
                    +2*ell*a*tt*k*t*mp.diff(hp, t)*nu-2*(hp(radius(t))-hp(t)) \
                    +2*a*t**2*tt*br*lt+2*a*t**2*tt*tau*ell**2*k*gv2+2*a*t**2*gt2
                norm = amp**2*min(t-1, 1)
                target = (product-ell*star)/norm
                erra = abs((areal-product)/norm)
                errd = abs(decomposition/norm-target)
                label = f"h{hi}_t{ti}_ell{ei}"
                check("areal_product_"+label, erra < mp.mpf('1e-30')*(1+abs(product/norm)), mp.nstr(erra, 8))
                check("finite_remainder_"+label, errd < mp.mpf('1e-30')*(1+abs(target)), mp.nstr(errd, 8))
                if ti == 2 and ei == 2:
                    check("positive_gradient_term_detected_h"+str(hi), abs(2*a*t**2*gt2/norm) > mp.mpf('1e-16'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 matched residual and inverse composition",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        exact_orbit_algebra()
        source_and_regional_bounds()
        forcing_norm_and_composition()
        weighted_projection_fixture()
        manufactured_orbit_fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S166 matched residual and inverse composition",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "mpmath": mp.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_residual_physical_forcing_bound": good,
        "compatible_residual_inverse_composition": good,
        "S166_linear_estimate_complete": good,
        "candidate_Hessian_inverse": False, "signed_compatibility_evaluated": False,
        "useful_Newton_certificate": False, "goal_complete": False,
        "unresolved": ["useful projected correction and nonlinear family/remainder",
                       "signed compatibility", "physical requirements"],
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
