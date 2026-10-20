"""S165: exact geometric identities and budgets; full metric inverse remains open."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "tt-form-and-scalar-fourth-graph"
OWN = ["checks/tt_form_and_scalar_fourth_graph.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-curved-quantum-tt", 157),
    ("stationary-weighted-tt-inverse", 52),
    ("stationary-scalar-outgoing-inverse", 107),
    ("global-compatible-scalar-response", 344),
]
PROTOCOL = {
    "arithmetic": "exact rational and symbolic; no fitted constants",
    "geometry": "independent 4D Christoffel and rough Laplacian calculation",
    "bolt_fixture": "a=rho, r=constant, polynomial c with exact regular b and eta",
    "pi_bracket": ["3", "22/7"],
    "parameter_box": {"x_min": 1000, "delta": [0, 1], "abs_d_max": 1},
    "benchmark": {"abs_d_over_x": "1.3e-22", "g_upper": "2000^-10",
                  "magnetic_v_upper": "4e-27/180"},
    "fixture_is_physical_preparation": False,
}
ROWS = []
R = S.Rational


def check(name, ok, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    result = S.factor(S.expand(S.simplify(expression)))
    check(name, result == 0, "" if result == 0 else result)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return data.get("runs", [data])[-1]


def hashes_of(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def input_paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
        out.update(hashes_of(latest(stem)))
    return sorted(out)


def predecessors():
    for stem, expected in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed"
              and len(run["checks"]) == expected and all(x["passed"] for x in run["checks"])
              and hashes_of(run) == hashes_of(run, True))
        for name, digest in hashes_of(run).items():
            check(stem+"_input_"+name, sha256((ROOT/name).read_bytes()).hexdigest() == digest)


def geometric_identities():
    rho = S.symbols("rho", positive=True)
    A, RR, b, c, eta, c1, c2, F2, F2p = S.symbols("A R b c eta c1 c2 F2 F2p")
    rules = {A: F2/2-A*A, RR: A*RR-RR*RR,
             b: -2*A*(b+c)-2*RR*(b-c), c: c1, c1: c2,
             eta: A*eta+b+c, F2: F2p}
    diff = lambda p: S.expand(sum(S.diff(p, v)*d for v, d in rules.items()))
    K = F2/2-2*A*RR+RR*RR
    T = -(b+3*c)/2-(A-RR)*eta
    first = -3*c1/2+2*RR*b-K*eta
    second = -3*(c2+A*c1)/2-(6*RR*RR+K)*b+(4*RR*RR-4*A*RR-K)*c-(diff(K)+2*A*K)*eta
    zero("return_first_identity", diff(T)-first)
    zero("return_second_identity", diff(diff(T))+A*diff(T)-second)
    lap_b = diff(rules[b])+(A+2*RR)*rules[b]-4*A*A*(b+c)-4*RR*RR*(b-c)
    zero("TT_radial_rough_laplacian_substitution",
         lap_b-(-F2*(b+c)-2*RR*RR*(b-c)-2*(A-RR)*c1))
    phi, hh = S.Function("phi")(rho), S.Function("H")(rho)
    theta = S.Function("theta")(rho)
    lap = lambda p: S.diff(p, rho, 2)+hh*S.diff(p, rho)
    chi = [S.cos(theta), S.sin(theta)]
    zero("unlocalization_identity",
         sum(t*lap(t*phi) for t in chi)+S.diff(theta, rho)**2*phi-lap(phi))
    radius, cc, beta, t2 = S.symbols("r C beta T2", positive=True)
    zero("orbit_margin_physical_measure",
         cc*(beta*cc/radius**2)*(radius**2/cc*t2)**2-beta*radius**2*t2**2)
    tt, E = S.symbols("T E")
    zero("fourth_graph_equation", (tt+4)*(tt-E)+(E-4)*tt+4*E-tt**2)

    # Direct coordinate calculation, independent of the frame formulas.
    time, angle, azimuth = S.symbols("time theta azimuth", real=True)
    coords = [time, rho, angle, azimuth]
    aa, rr, bb, cc = [S.Function(n)(rho) for n in ["a", "r", "b", "c"]]
    metric = [aa**2, S.Integer(1), rr**2, rr**2*S.sin(angle)**2]
    inv = [1/q for q in metric]
    frame = [-bb-2*cc, bb, cc, cc]
    H = [metric[i]*frame[i] for i in range(4)]
    Gamma = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                expr = ((S.diff(metric[j], coords[i]) if k == j else 0)
                        +(S.diff(metric[i], coords[j]) if k == i else 0)
                        -(S.diff(metric[i], coords[k]) if i == j else 0))*inv[k]/2
                Gamma[k, i, j] = S.simplify(expr)
    nH = {}
    for j in range(4):
        for a in range(4):
            for d in range(4):
                nH[j, a, d] = S.simplify(
                    (S.diff(H[a], coords[j]) if a == d else 0)
                    -Gamma[d, j, a]*H[d]-Gamma[a, j, d]*H[a])
    Avec, Rvec = S.diff(aa, rho)/aa, S.diff(rr, rho)/rr
    scal = lambda p: S.diff(p, rho, 2)+(Avec+2*Rvec)*S.diff(p, rho)
    target_b = scal(bb)-4*Avec**2*(bb+cc)-4*Rvec**2*(bb-cc)
    target_c = scal(cc)+2*Rvec**2*(bb-cc)
    targets = [-target_b-2*target_c, target_b, target_c, target_c]
    for a in range(4):
        for d in range(4):
            value = 0
            for i in range(4):
                v = S.diff(nH[i, a, d], coords[i])
                v -= sum(Gamma[k, i, i]*nH[k, a, d]
                         +Gamma[k, i, a]*nH[i, k, d]
                         +Gamma[k, i, d]*nH[i, a, k] for k in range(4))
                value += inv[i]*v
            expected = targets[a]*metric[a] if a == d else 0
            zero(f"coordinate_rough_laplacian_{a}_{d}", value-expected)
    zero("rough_laplacian_trace", sum(targets))
    scalar_coord = sum(inv[i]*(S.diff(phi, coords[i], 2)
                       -sum(Gamma[k, i, i]*S.diff(phi, coords[k]) for k in range(4)))
                       for i in range(4))
    zero("coordinate_scalar_laplacian", scalar_coord-scal(phi))


def polynomial_identities_and_envelopes():
    u, t, h, delta, d = S.symbols("u t h delta d", real=True)
    # h=1/x; all polynomial majorants use |h|<=1/1000, |d|,|delta|<=1.
    f = 1+h*u**2-(4+(2+delta)*h)*u**3+(3+(1+d)*h)*u**4-2*d*h*u**5+d*h*u**6
    z = d*h*(4*u**5-6*u**6)
    j = -(3+(1+d)*h)*u**4+5*d*h*u**5-R(9, 2)*d*h*u**6
    w = -(2+(1+delta/2)*h)*u**3+(3+(1+d)*h)*u**4-R(10, 3)*d*h*u**5+R(5, 2)*d*h*u**6
    Lop = lambda p: S.expand(f*u**2*S.diff(p, u, 2)+(u**2*S.diff(f, u)-2*f*u)*S.diff(p, u))
    Ft = t**2*f.subs(u, 1/t)  # F/x, t=r/B.
    K = S.expand(S.diff(Ft, t, 2)/2-S.diff(Ft, t)/t+Ft/t**2).subs(t, 1/u)
    Kpoly = h*u**2-6*(2+(1+delta/2)*h)*u**3+6*(3+(1+d)*h)*u**4-20*d*h*u**5+15*d*h*u**6
    zero("seed_K_polynomial", K-Kpoly)
    v = S.Function("v")(u)
    vv = v.subs(u, 1/t)
    direct = S.diff(Ft*S.diff(vv, t), t)+2*Ft*S.diff(vv, t)/t
    zero("seed_scalar_radial_operator", direct.subs(t, 1/u)-Lop(v))
    e, mag = S.symbols("e v")
    Ebar = (1+11*e/3)*z+e*Lop(z)-11*e*z*z/72+28*e*w*w+R(22, 3)*e*j*j-10*mag*u**4
    raw = -12+z+12*(1+R(11, 6)*e)+e*(Lop(z)-R(11, 72)*(-12+z)**2
                                                   +R(7, 12)*48*w*w+R(11, 6)*4*j*j)-10*mag*u**4
    zero("complete_local_trace_polynomial", raw-Ebar)
    check("local_trace_minimum_decay", min(mon[0] for mon, coeff in S.Poly(S.expand(Ebar), u).terms()) >= 4)
    check("local_trace_maximum_degree", S.degree(Ebar, u) <= 12)

    def coeff_majorant(p):
        terms = S.Poly(S.expand(p), u, h, delta, d).terms()
        return sum(abs(coef)*R(1, 1000)**mon[1] for mon, coef in terms)
    for name, p, bound in [
        ("f", f, 9), ("ufprime", u*S.diff(f, u), 50),
        ("u2fsecond", u*u*S.diff(f, u, 2), 250),
        ("j", j, R("3.6")), ("ujprime", u*S.diff(j, u), 20),
        ("w", w, R("5.1")), ("uwprime", u*S.diff(w, u), 40),
        ("K", Kpoly, 31), ("uKprime", u*S.diff(Kpoly, u), 109),
        ("Box_j", Lop(j), 2000),
    ]:
        got = coeff_majorant(p)
        check("polynomial_envelope_"+name, got < bound, got)
    for name, p, bound in [
        ("z", z/(d*h), 10), ("uzprime", u*S.diff(z,u)/(d*h), 56),
        ("u2zsecond", u*u*S.diff(z,u,2)/(d*h), 260),
        ("Box_z", Lop(z)/(d*h), 7000),
    ]:
        got = coeff_majorant(S.cancel(p))
        check("polynomial_envelope_"+name, got <= bound, got)
    check("K_proper_derivative", 5*109**2 < 300**2)
    check("Box_z_degree", S.degree(Lop(z), u) <= 12)
    check("trace_zero_Einstein_limit", S.simplify(Ebar.subs({h:0, d:0, mag:0}).subs(u,0)) == 0)


def bolt_diagnostic():
    rho = S.symbols("rho", positive=True)
    c0, c2, c4, eta1 = S.symbols("c0 c2 c4 eta1")
    c = c0+c2*rho**2+c4*rho**4
    b = -c0-c2*rho**2/2-c4*rho**4/3
    eta = eta1*rho+c2*rho**3/4+c4*rho**5/6
    A = 1/rho
    T = -(b+3*c)/2-A*eta
    lap = lambda p: S.diff(p,rho,2)+A*S.diff(p,rho)
    zero("bolt_TT_equation", S.diff(b,rho)+2*A*(b+c))
    zero("bolt_return_equation", S.diff(eta,rho)-A*eta-b-c)
    zero("bolt_return_first", S.diff(T,rho)+R(3,2)*S.diff(c,rho))
    zero("bolt_return_second", lap(T)+R(3,2)*lap(c))
    zero("bolt_rough_radial", lap(b)-4*A*A*(b+c)+2*A*S.diff(c,rho))
    zero("bolt_scalar_graph", lap(c)-(4*c2+16*c4*rho**2))
    zero("bolt_shape_graph", lap(T)+6*c2+24*c4*rho**2)
    zero("bolt_A_cprime_limit", S.limit(A*S.diff(c,rho),rho,0)-2*c2)
    zero("bolt_A2_d1_limit", S.limit(A*A*(b+c),rho,0)-c2/2)
    zero("bolt_A_eta_limit", S.limit(A*eta,rho,0)-eta1)


def budgets():
    pi = R(22,7)
    alpha = S.Integer(2000)**-5
    d_x, emag, ev = R("1.3e-22"), alpha**2/180, R("4e-27")/180
    check("unlocalization_loss", R(500**4,300)*R("3.4e-37") < R("1e-28"))
    check("retained_energy_margin", 1-R(".004")-R("1e-28")-R("1e-6") > R(".99"))
    check("orbit_graph_constant", 800*pi*pi < 89**2)
    for name, val, bound in [
        ("c_L2", 8*pi/R(".99"), R("5.1")**2),
        ("c_gradient", 320*pi/R(".99"), 34**2),
        ("h_L2", 176*pi/R(".99"), 24**2),
        ("h_gradient", 41600*pi/R(".99"), 364**2),
        ("eta_L2", 144*pi/R(".99"), 22**2),
    ]:
        check(name+"_constant", val < bound)
    check("return_divergence_factor", 4*90 < 19**2)
    check("return_divergence_constant", 28*R("5.1") < 150)
    check("scalar_graph_lower_terms", R(2,3)*(244+300*R("4.25")+2*40*3) < 1300)
    check("return_eta_factor", 18 < R("4.25")**2)
    check("scalar_Bochner_reserve", 14 < 4**2)
    check("tensor_laplacian_value_reserve", 2*(21*3+2*5*3)+3*(2*5*3) == 276)
    check("tensor_laplacian_gradient_reserve", 2*(2*4+R("4.5")) < 26)
    check("tensor_laplacian_gradient_root", 4*5 < R("4.5")**2)
    check("curvature_derivative_bound", 1000+4*100+R(2,3) < 2000)
    check("tensor_Bochner_gradient", 14+4*54+4000 < 10000)
    check("tensor_Bochner_value", 4000 < 10000)
    check("scalar_graph_constant", R(2,3)*89+(1300*R("5.1")+5*34)*alpha < 61)
    check("rough_tensor_graph_constant", 7*61+(26*34+276*R("5.1"))*alpha < 430)
    check("tensor_Hessian_constant", 430+(100*364+100*24)*alpha < 500)
    check("Box_j_generic_envelope", (18*9+6*50)*R("3.6") < 2000)
    check("Box_z_generic_envelope", 9*260+(50+18)*56 < 7000)
    check("trace_gradient_polynomial_reserve",
          R("2.25")*(56*(1+11*emag/3)*d_x+84000*emag*d_x
             +R(11,36)*emag*10*56*d_x**2
             +28*emag*2*R("5.1")*40+R(22,3)*emag*2*R("3.6")*20+40*ev)
          < 200*d_x+100000*emag+100*ev)
    check("trace_gradient_bound", 200*d_x+100000*emag+100*ev < R("4e-20"))
    check("abar_gradient_bound", R("2.25")*56*d_x/36+R(11,9)*30 < 40)
    check("mixed_second_derivative_coefficient", R(7,12)*6*40 == 140)
    check("mixed_rough_derivative_coefficient", R(11,6)*7 < 13)
    algebra = R(7*11*13,36)+4+R(7,12)*12320+R(11,6)*(7*(2*54+2*14)+R(13,2)*7+2*7**2)
    check("mixed_algebraic_reserve", algebra < 20000, algebra)
    # Divide the vertex inequality by 1/sqrt(g) before testing its constant.
    qcoef = 14*61+140*500+13*430
    small = 240*34+8000*R("5.1")+20000*24
    vertex = R(9,8)*(7*24+alpha/180*(qcoef+small*alpha)
                     +20*ev*24+R("6e-22")*150+R("4e-20")*22)
    check("strong_mixed_vertex_constant", vertex < 250, vertex)
    check("auxiliary_source_bound", 3*30000+6100*R(".2")*alpha < 91000)
    check("strong_Schur_source_constant", 80*pi**2*250**2*(8*pi/R(".99")) < 36000**2)
    check("strong_trace_source_constant", 80*pi**2 < 800)
    # g^(3/4)<alpha*(6e-9), using S164's exact root bound.
    check("quarter_root_g", alpha < R("6e-9")**2)
    local = 5*R("3.2e9")*alpha+40*R("4e8")*alpha*R("6e-9")+40*pi*d_x*R("4.4e7")
    check("local_scalar_strong_source", local < R("1e-6"), local)
    check("total_scalar_source", 36000+250*91000+R("1e-6") < R("2.4e7"))
    check("E_upper_bound", 60*pi+R(11,3)*alpha**2 < 190)
    check("fourth_scalar_graph_budget",
          R("2.4e7")*alpha+190*R("3.2e9")+4*190*R("4.4e7")*alpha < R("6.2e11"))
    check("full_shape_form_bound", (3500-R("5.1"))**2*80*3**2 > 91000**2)
    check("full_shape_orbit_graph_bound", 89*3500 < R("3.2e5"))
    check("full_shape_Hessian_bound", 500*3500 < R("1.8e6"))


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = input_paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before, exception = hashes(), None
    print(json.dumps({"event":"start", "story":"S165", "input_count":len(ps)}), flush=True)
    try:
        predecessors()
        geometric_identities()
        polynomial_identities_and_envelopes()
        bolt_diagnostic()
        budgets()
    except Exception:
        exception = traceback.format_exc()
        print(exception, file=sys.stderr, flush=True)
    after = hashes()
    passed = exception is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story":"S165", "created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if passed else "failed", "python":sys.version, "sympy":S.__version__,
        "protocol":PROTOCOL, "input_sha256_before":before, "input_sha256_after":after,
        "unchanged_inputs":before == after, "checks":ROWS, "exception":exception,
        "scope":"actual TT form regularity, strong mixed source and compatible scalar fourth graph",
        "unresolved":["fast fourth/log mass-extended shape inverse","source/reconstruction constants",
                      "nonlinear family","actual signed preparation source","physical completion"],
    }
    if not passed:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs":[]}
    data["runs"].append(run)
    dest.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"], "checks":len(ROWS), "unchanged_inputs":before == after}),flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
