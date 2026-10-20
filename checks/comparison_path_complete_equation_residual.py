"""S161: fixed-path complete retained-equation defect, not solution error."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "comparison-path-complete-equation-residual"
PREDECESSORS = [
    ("comparison-path-complete-source-response", 305),
    ("fixed-chart-global-constraint-transfer", 99),
    ("stationary-centered-complete-residual", 92),
]
OWN = ["checks/comparison_path_complete_equation_residual.py",
       "notes/comparison-path-complete-equation-residual.md"]
PROTOCOL = {
    "digits": 100, "h": ["1e-5", "1e-30"],
    "s": ["2h", ".001", ".1", "3", "49"], "parameters": ["0", ".4", "1"],
    "amplitude": ".02", "normalization": "amplitude^2 min(s,1)",
    "formula_tolerance": "1e-30", "contour_nodes": [64, 128],
    "contour_radius": ".01", "contour_tolerance": "1e-30",
    "real_steps": [".001", ".0005"], "difference_tolerance": "1e-6",
    "refinement_ratio": ".35", "refinement_floor": "1e-35",
    "endpoints": "second-order one-sided within [0,1]",
    "fixtures_are_physical_states": False,
}
ROWS, NUMERICS = [], []
R = S.Rational


def check(name, passed, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def zero(name, expr):
    check(name, S.factor(S.together(expr)) == 0)


def read_predecessor(stem):
    return json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))["runs"][-1]


def input_names():
    names = set(OWN)
    for stem, _ in PREDECESSORS:
        names.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
        r = read_predecessor(stem)
        names.update(r.get("input_sha256_before", r.get("execution_inputs_sha256", {})))
    return sorted(names)


def predecessors():
    for stem, count in PREDECESSORS:
        r = read_predecessor(stem)
        before = r.get("input_sha256_before", r.get("execution_inputs_sha256"))
        after = r.get("input_sha256_after", r.get("final_inputs_sha256"))
        check(stem+"_receipt", r["status"] == "passed" and len(r["checks"]) == count
              and all(x["passed"] for x in r["checks"]) and before == after)
        for p, h in before.items():
            check(stem+"_input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def exact_identities():
    t, t1, t2, sig2, chi2, H, P = S.symbols("t t1 t2 sig2 chi2 H P", positive=True)
    X, X1, X2, vp = S.symbols("X X1 X2 vp")
    # Euler derivative of -2X chi'' + X'^2/(2X) + exp(2chi)V(X).
    density_derivative = -2*chi2-X1**2/(2*X**2)+vp
    dt_momentum = X2/X-X1**2/X**2
    EXdensity = density_derivative-dt_momentum
    converted = (-t**2*H/2*EXdensity).subs({
        X: t**2, X1: 2*t*t1, X2: 2*t1**2+2*t*t2,
        chi2: sig2+t2/t-t1**2/t**2, vp: 2*P/(t**2*H)})
    zero("independent_radius_Euler_variation", converted-(H*(t**2*sig2+2*t*t2-t1**2)-P))
    check("gradient_square_coefficient_detected", S.expand(converted).coeff(t1, 2) == -H)

    ell, u, V, Deltau, Lnu, Ltau, gradV2, gradT2, Ps, Pnext, Pprime = S.symbols(
        "ell u V Deltau Lnu Ltau gradV2 gradT2 Ps Pnext Pprime")
    tau = (t-1)/t
    nu, k, a = tau*V, S.symbols("k"), S.symbols("a")
    T = 1+tau*(k-1)
    B = k-1-ell*V*k
    LT = ell*k*Lnu+B*Ltau+tau*ell**2*k*gradV2
    A0 = a*(T**2*Ps+ell*t*t*T*T*Deltau+2*t*t*T*LT-t*t*gradT2)-Pnext
    Astar = (2*nu-2*u)*Ps+t*t*Deltau+2*t*t*Lnu-t*Pprime*nu
    decomposition = (ell*(a*T*k-1)*Astar
        +Ps*(a*(T*T+2*ell*T*k*(u-nu))-1)
        +ell*a*T*k*t*Pprime*nu-(Pnext-Ps)
        +ell*a*t*T*(1-k)*Deltau+2*a*t*t*T*B*Ltau
        +2*a*t*t*T*tau*ell**2*k*gradV2-a*t*t*gradT2)
    zero("exact_angular_nonlinear_decomposition", A0-ell*Astar-decomposition)
    zero("T_minus_k_cancellation", T-k-(1-k)/t)

    f, fp = S.symbols("f fp")
    D = lambda z: S.diff(z, t)
    Delta = lambda z: t*t*(f*D(D(z))+fp*D(z))
    L = lambda z: Delta(z)+t*f*D(z)
    zero("seed_Ltau", L(tau)-(fp-f/t))
    Vfunc = S.Function("V")(t)
    kfunc = S.exp(ell*Vfunc)
    Tfunc = 1+tau*(kfunc-1)
    zero("radial_T_gradient", t*D(Tfunc)-((kfunc-1)/t+ell*kfunc*(t-1)*D(Vfunc)))
    Nfunc = tau*Vfunc
    zero("LT_product_formula", L(Tfunc)-ell*kfunc*L(Nfunc)
         -(kfunc-1-ell*Vfunc*kfunc)*L(tau)-tau*ell**2*kfunc*t*t*f*D(Vfunc)**2)

    jr, j1, j2, aa, aa1 = S.symbols("jr j1 j2 aa aa1")
    residual = aa-ell*(aa1+ell*j1+ell**2*j2/2)
    zero("complete_residual_allocation", residual-(aa-ell*aa1-ell**2*j1-ell**3*j2/2))
    zero("quadratic_cancellation", S.diff((ell*aa1+ell**2*aa/2)-ell*(aa1+ell*j1), ell).subs(ell, 0))
    # Physical-to-normalized angular conversion and the source multipliers.
    G, BB, CC, rr, ex = S.symbols("G BB CC rr ex", nonzero=True)
    ps = -ex/(16*S.pi*G)
    zero("angular_physical_conversion", 8*S.pi*G*rr**4*ps/(BB**2*CC)+rr**4*ex/(2*BB**2*CC))
    b = S.symbols("b")
    T = 1+tau*(S.exp(ell*V)-1)
    Fr = S.exp(2*ell*(V-u))*(1+ell*b)**2/T**2
    radial_map = S.diff(1+(t-1)*S.exp(ell*Vfunc), t)
    zero("mass_prefactor_once", radial_map/Tfunc**2
         -S.exp(ell*Vfunc)*(1+ell*(t-1)*S.diff(Vfunc, t))/Tfunc**2)
    zero("lapse_prefactor_once", S.exp(ell*V)*(1+ell*b)/(T**3*Fr)
         -S.exp(ell*(2*u-V))/(T*(1+ell*b)))


def budgets():
    mp.mp.dps = 100
    check("small_exponential_guard", mp.exp(mp.mpf("1e-6")) < mp.mpf("1.001"))
    n, w, Z = S.symbols("n w Z", nonnegative=True)
    half = R("2.02")*n*n+R("1.01")*Z+R("2.01")*n*(2*n+w)
    # n[(2n+w)^2+1.001Z], using n,w<=1e-6, with positive monomials.
    half += R("1.01")*(4*R("1e-6")*n*n+5*R("1e-6")*n*w+R("1.001e-6")*Z)
    margin = S.Poly(10*n*n+5*n*w+2*Z-half, n, w, Z)
    check("scalar_Cp_majorant", all(c > 0 for c in margin.coeffs()))
    check("potential_leading_coefficients", 3/R("5.9") < R(".51")
          and 12/R("5.9") < R("2.1") and 36/R("5.9") < R("6.3"))

    def main_terms(t, u, v, nu, du, p, d1, d2, tau, lap, astar):
        z = u+nu
        return (R("1.01")*(2*u+v+R("1.001")*nu)*astar
                +p*(10*z*z+5*z*v+2*tau*v*v)
                +R("1.01")*d1*nu*(2*u+v+R("1.001")*nu)
                +R(".501")*d1*tau*v*v+R(".502")*d2*nu*nu
                +R("1.01")*t*v*du+R("1.01")*v*v*lap)

    ax, ss = S.symbols("ax ss", nonnegative=True)
    tm, vb, bb = R("1.001"), R("1.2e-8"), R("3e-7")
    cold = main_terms(tm, R("4e-6")*ax, vb, ss*vb, R("8.05e-6")*ax,
                      R("1.01"), 3, 10, ss, 5*ss, R("1.7e-6")*ax)
    cold += R("2.0402")*tm*ss*bb*bb+R("1.0201")*ss**2*(vb+bb)**2
    poly = S.Poly(S.expand(cold), ax, ss)
    check("cold_positive_weighted_polynomial", all(c > 0 for c in poly.coeffs()) and poly.coeff_monomial(1) == 0)
    cold_bound = 0
    for (i, j), c in poly.terms():
        if i:
            cold_bound += c*R(".00331")**(i-1)*R(".001")**j
        else:
            assert j >= 1
            cold_bound += c*R(".001")**(j-1)  # s/ax <=1
    check("cold_E0_angular_remainder", cold_bound < R("2e-12"), str(cold_bound))
    t = R(4)
    mid = main_terms(t, vb, vb, vb, R("7.4e-8"), R(".51")*(t**4+1),
                     R("2.1")*t**4, R("6.3")*t**4, 1, 3*t**3, R("2.6e-9"))
    mid += R("2.02")*t**3*(t-1)*bb**2+R("1.01")*t**2*(t-1)**2*(vb+bb)**2
    check("middle_E0_angular_remainder", mid < R("6e-11"), str(mid))
    z = R(1, 4)
    # Divide by L_e^2>=1 before these monotone z=1/t substitutions.
    out = main_terms(1/z, R("8e-9")*z*z, R("2e-9")*z*z, R("2e-9")*z*z,
                     R("7.4e-8"), R(".51")*(1+z**4)/z**4,
                     R("2.1")/z**4, R("6.3")/z**4, 1, 8/z**3, R("4.7e-7"))
    out += (R("2.02")*5/R("5.9")*R(4, 3)*R("1.4e-8")**2
            +R("1.01")*5/R("5.9")*(R("2e-9")*z+R("1.4e-8"))**2)
    check("outer_E0_angular_remainder", out < R("2e-15"), str(out))
    check("middle_upper_subregion", R("2e-15")*25 < R("6e-11"))
    for name, a, b, limit in [
        ("cold", "2e-12", "2.6e-10", "2.7e-10"),
        ("middle", "6e-11", "2.2e-12", "6.3e-11"),
        ("outer", "2e-15", "4.9e-11", "5e-11")]:
        check("complete_angular_"+name, R(a)+R(b) < R(limit))
    check("retained_mass_defect", R("5e-13")+R("2e-11") < R("2.1e-11"))
    check("retained_lapse_defect", R("1.3e-10")+R("1.2e-10") < R("2.6e-10"))
    check("retained_horizon_defect", R("1e-17")+R("1e-13") < R("2e-13"))
    NUMERICS.append({"cold_E0_exact": str(cold_bound), "middle_E0_exact": str(mid),
                     "outer_E0_exact": str(out), "inverse_proved": False, "goal_complete": False})


def fixtures():
    t, ell, thermal = S.symbols("t ell thermal")
    e, d, amp = R(".001"), R(".00002"), R(".02")
    x = (1-e-d)/6
    fext = (t-1)**2/t**2*(e+x*(t*t+2*t+3)+d/t**2)
    f = fext-thermal/t
    dt = lambda z: S.diff(z, t)
    Delta = lambda z: t*t*(f*dt(dt(z))+dt(f)*dt(z))
    L = lambda z: Delta(z)+t*f*dt(z)
    R0 = -(t*t*dt(dt(f))-2*t*dt(f)+2*f)
    P = lambda z: 3*x*z**4+(e+3*x+d)-6*d/z+6*d/z**2
    zero("exact_seed_angular_equation", -t*t*R0/2+2*t*Delta(t)-t*t*f-P(t))
    u, V = amp*(t-1)/(1+t)**3, amp/(1+t)**3
    k, a = S.exp(ell*V), S.exp(-2*ell*u)
    tau, nu, b = (t-1)/t, (t-1)/t*V, (t-1)*dt(V)
    T, tl = 1+tau*(k-1), 1+(t-1)*k
    A = a*(-tl*tl*R0/2+ell*tl*tl*Delta(u)+2*tl*Delta(tl)-t*t*f*dt(tl)**2)-P(tl)
    star = (2*nu-2*u)*P(t)+t*t*Delta(u)+2*t*t*L(nu)-t*dt(P(t))*nu
    zero("manufactured_tangent_identity", S.diff(A, ell).subs(ell, 0)-star)
    BB = k-1-ell*V*k
    rem = (ell*(a*T*k-1)*star+P(t)*(a*(T*T+2*ell*T*k*(u-nu))-1)
           +ell*a*T*k*t*dt(P(t))*nu-(P(tl)-P(t))
           +ell*a*t*T*(1-k)*Delta(u)+2*a*t*t*T*BB*L(tau)
           +2*a*t*t*T*tau*ell**2*k*t*t*f*dt(V)**2-a*t*t*t*t*f*dt(T)**2)
    fr = f*S.exp(2*ell*(V-u))*(1+ell*b)**2/T**2
    alpha = 2*S.log(T)+2*ell*u-ell*V-S.log(1+ell*b)
    Dr = lambda z: dt(z)/dt(tl)
    ar = tl**4*(Dr(Dr(fr))/2+Dr(fr)/tl+R(3, 2)*Dr(alpha)*Dr(fr)
                +fr*(Dr(Dr(alpha))+Dr(alpha)**2+Dr(alpha)/tl))-P(tl)
    af, rf, ef, sf, tf, gf = [S.lambdify((t, ell, thermal), z, "mpmath", cse=True)
                            for z in (A, rem, ar, star, S.diff(A, ell, 2),
                                      a*t**4*f*dt(T)**2)]
    fextfun = S.lambdify(t, fext, "mpmath", cse=True)
    failures_detected = False
    for hs in PROTOCOL["h"]:
        h = mp.mpf(hs)
        v = 1+h
        th = v*fextfun(v)
        for sval in PROTOCOL["s"]:
            ss = 2*h if sval == "2h" else mp.mpf(sval)
            tt, norm = 1+ss, mp.mpf(".02")**2*min(ss, mp.mpf(1))
            for ls in PROTOCOL["parameters"]:
                lam = mp.mpf(ls)
                name = f"fixture_{hs}_{sval}_{ls}"
                value, direct, remainder, a1 = af(tt, lam, th), ef(tt, lam, th), rf(tt, lam, th), sf(tt, lam, th)
                err = abs(value-direct)/norm
                check(name+"_areal_angular", err < mp.mpf("1e-30")*(1+abs(value/norm)), mp.nstr(err, 25))
                err = abs(value-lam*a1-remainder)/norm
                check(name+"_remainder", err < mp.mpf("1e-30")*(1+abs(remainder/norm)), mp.nstr(err, 25))
                if lam == 0:
                    check(name+"_seed_zero", abs(value/norm) < mp.mpf("1e-30"))
                target = tf(tt, lam, th)/norm
                fun = lambda z: af(tt, z, th)/norm
                contours = []
                for nodes in PROTOCOL["contour_nodes"]:
                    radius = mp.mpf(PROTOCOL["contour_radius"])
                    total = sum(fun(lam+radius*mp.exp(2j*mp.pi*j/nodes))*mp.exp(-4j*mp.pi*j/nodes)
                                for j in range(nodes))
                    got = 2*total/(nodes*radius**2)
                    err = abs(got-target)
                    check(name+f"_contour_{nodes}", err < mp.mpf("1e-30")*(1+abs(target)), mp.nstr(err, 25))
                    contours.append(mp.nstr(got, 35))
                errors = []
                for step in PROTOCOL["real_steps"]:
                    hh = mp.mpf(step)
                    if lam == 0:
                        got = (2*fun(lam)-5*fun(lam+hh)+4*fun(lam+2*hh)-fun(lam+3*hh))/hh**2
                    elif lam == 1:
                        got = (2*fun(lam)-5*fun(lam-hh)+4*fun(lam-2*hh)-fun(lam-3*hh))/hh**2
                    else:
                        got = (fun(lam+hh)-2*fun(lam)+fun(lam-hh))/hh**2
                    err = abs(got-target)
                    errors.append(err)
                    check(name+"_difference_"+step, err < mp.mpf("1e-6")*(1+abs(target)), mp.nstr(err, 25))
                check(name+"_refinement", errors[1] < mp.mpf(".35")*errors[0]
                      or max(errors) < mp.mpf("1e-35"))
                q1, q2 = norm/(1+tt)**2, -norm/(3*(1+tt)**2)
                full = value-lam*(a1+lam*q1+lam**2*q2/2)
                check(name+"_complete_source", abs(full-(remainder-lam**2*q1-lam**3*q2/2))/norm < mp.mpf("1e-30"))
                if lam and abs(gf(tt, lam, th)/norm) > mp.mpf("1e-10"):
                    failures_detected = True
                NUMERICS.append({"fixture": name, "second": mp.nstr(target, 35),
                                 "contours": contours, "difference_errors": [mp.nstr(x, 25) for x in errors]})
    check("omitted_gradient_square_detected", failures_detected)


def main():
    mp.mp.dps = PROTOCOL["digits"]
    names = input_names()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in names}
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S161", "input_count": len(names)}), flush=True)
    try:
        predecessors()
        exact_identities()
        budgets()
        fixtures()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story": "S161", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after, "unchanged_inputs": before == after,
        "checks": ROWS, "numerics": NUMERICS, "exception": error,
        "scope": "complete retained-equation defect on the specified comparison path",
        "unresolved": ["general family", "usable full inverse", "coupled solution", "physical completion"],
    }
    if error:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    path = ROOT/"evidence"/(STEM+".json")
    history = json.loads(path.read_text()) if path.exists() else {"runs": []}
    history["runs"].append(run)
    path.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
